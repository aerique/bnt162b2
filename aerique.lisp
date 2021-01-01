;;;; aerique.lisp
;;;;
;;;; [My](https://github.com/aerique/bnt162b2) attempts at finding an
;;;; algorithm to get from the virus RNA to the vaccine RNA.
;;;;
;;;; Requires Common Lisp ([SBCL](http://www.sbcl.org/) is available for most
;;;; Unix's) and [Quicklisp](https://www.quicklisp.org/).
;;;;
;;;; Assuming you've never used Common Lisp but managed to install `sbcl`:
;;;;
;;;; - follow the steps at https://www.quicklisp.org/beta/#installation
;;;  - OR:
;;;;     1. `sbcl --load quicklisp.lisp --eval "(progn (quicklisp-quickstart:install) (quit))"`
;;;;     2. `echo "(load (merge-pathnames \"quicklisp/setup.lisp\" (user-homedir-pathname)))" > $HOME/.sbclrc`
;;;;
;;;; Usage: `sbcl --load aerique.lisp --eval "(main)"`

;;; Packages

(ql:quickload :cl-csv :silent t)


;;; Common Functions

(defun mk-pathname (pathname-or-string)
  (if (typep pathname-or-string 'string)
      (merge-pathnames pathname-or-string)
      pathname-or-string))


(defun print-hash-table (hash-table)
  (maphash (lambda (k v)
             (format t "~S: ~S~%" k v))
           hash-table))


;;; CSV Functions

(defun print-side-by-side (parsed-side-by-side)
  "PARSED-SIDE-BY-SIDE is the output of READ-SIDE-BY-SIDE-CSV."
  (format t "abspos | codonOrig | codonVaccine~%~
             -------+-----------+-------------~%")
  (loop for codon across parsed-side-by-side
        do (format t "  ~4D |       ~A | ~A~%"
                   (elt codon 0) (elt codon 1) (elt codon 2))))


(defun read-side-by-side-csv (&optional (csv "side-by-side.csv"))
  (let ((pathname (if (typep csv 'string)
                      (merge-pathnames csv)
                      csv)))
    ;; parse into vector for speed
    (loop with result = (make-array 1273 :element-type 'vector :fill-pointer 0
                                         :initial-element (vector))
          for row in (rest (cl-csv:read-csv pathname))
          do (vector-push-extend ;; abspos codonOrig codonVaccine
                                 (vector (parse-integer (first row))
                                         (second row) (third row))
                                 result)
          finally (return result))))


(defun read-codon-table-csv (&optional (csv "codon-table-grouped.csv"))
  (let ((pathname (if (typep csv 'string)
                      (merge-pathnames csv)
                      csv))
        (ht (make-hash-table :test #'equal :size 64)))
    (loop for row in (rest (cl-csv:read-csv pathname))
          for aminoacid = (first row)
          for codon = (second row)
          do (setf (gethash codon ht) aminoacid))
    ht))


(defun make-equivalence-tables (&optional (csv "codon-table-grouped.csv"))
  "Returns a plist with two hashtables:
  - C2A: returns the aminoacid a codon is translated into,
  - A2C: returns the codons to create an aminoacid."
  (let ((c2a (make-hash-table :test #'equal :size 64))
        (a2c (make-hash-table :test #'equal :size 21)))
    (loop for row in (rest (cl-csv:read-csv (mk-pathname csv)))
          for aminoacid = (first row)
          for codon = (second row)
          do (setf (gethash codon c2a) aminoacid))
    (loop for row in (rest (cl-csv:read-csv (mk-pathname csv)))
          for aminoacid = (first row)
          for codon = (second row)
          do (setf (gethash aminoacid a2c)
                   (append (gethash aminoacid a2c) (list codon))))
    (list :c2a c2a :a2c a2c)))


;;; Functions

(defun print-strand-for-rnafold (codons)
  (loop for codon across codons
        append (list (elt codon 1) (elt codon 1) (elt codon 2)) into result
        finally (return (coerce result 'string))))


(defun compare-against-vaccine (codons parsed-side-by-side)
  "CODONS should be a vector of codons where each codon is a three letter
  string: #(\"TAA\" \"ACA\" ...)
  Returns the percentage of codons EQUAL to the codons in the vaccine."
  (loop with n_different = 0
        for codon across codons
        for i from 0
        do (when (string/= codon (elt (elt parsed-side-by-side i) 2))
             (incf n_different))
        finally (return (- 100 (* (/ n_different (length parsed-side-by-side))
                                  100.0)))))


(defun get-original-codons (parsed-side-by-side)
  "PARSED-SIDE-BY-SIDE is the output of READ-SIDE-BY-SIDE-CSV."
  (loop for row across parsed-side-by-side
        collect (elt row 1) into result
        finally (return (coerce result 'vector))))


;;; No Operation 'Algorithm

(defun nop (parsed-side-by-side &key verbose)
  (declare (ignore verbose))
  (get-original-codons parsed-side-by-side))


;;; Bert's Sample Algorithm: ../00-Resources/bnt162b2-git/3rd-gc.go

(defun berts-algorithm (parsed-side-by-side &key (verbose t))
  (let ((codons (get-original-codons parsed-side-by-side))
        (equiv (read-codon-table-csv)))
    (loop for codon across codons
          for i from 0
          for aminoacid = (gethash codon equiv)
          for new-codon = (copy-seq codon)
          do (unless (or (char= #\C (elt codon 2))
                         (char= #\G (elt codon 2)))
               (setf (elt new-codon 2) #\G)
               (if (string= aminoacid (gethash new-codon equiv))
                   (progn (when verbose
                            (format t "~4D: (->G) ~A = ~A (~A)~%"
                                    i new-codon codon aminoacid))
                          (setf (elt codons i) new-codon))
                   (progn (setf (elt new-codon 2) #\C)
                          (when (string= aminoacid (gethash new-codon equiv))
                            (when verbose
                              (format t "~4D: (->C) ~A = ~A (~A)~%"
                                      i new-codon codon aminoacid))
                            (setf (elt codons i) new-codon))))))
    codons))


;;; My Algorthim(s)

(defun n-gc (codon)
  "Returns the number of G and C in CODON."
  (loop for c across codon
        when (or (char= c #\C)
                 (char= c #\G))
          sum 1))


;; Attempt 1: replace a codon with a codon that has a higher G or C count,
;;            otherwise keep the current codon
;;
;; Specific: the codons belonging to an aminoacid are evaluated in the same
;;           order as they appear in `codon-table-grouped.csv`
;;
;; Result: 62.72%
(defun replace-with-higher-gc-codon-01 (parsed-side-by-side &key (verbose t))
  (let* ((codons (get-original-codons parsed-side-by-side))
         (tables (make-equivalence-tables))
         (c2a (getf tables :c2a))   ; codon to aminoacid
         (a2c (getf tables :a2c)))  ; aminoacid to codons
    (loop for codon across codons
          for i from 0
          for ac = (gethash codon c2a)      ; get aminoacid for this codon
          for ac-codons = (gethash ac a2c)  ; get all codons for the aminoacid
          for ngc = (n-gc codon)  ; number of Gs and Cs in `codon`
          ;; return first codon with highest number of Gs and Cs
          for highest = (loop with highest-codon = (first ac-codons)
                              with highest-gc = (n-gc highest-codon)
                              for ac-codon in (rest ac-codons)
                              for ac-codon-gc = (n-gc ac-codon)
                              do (when (> ac-codon-gc highest-gc)
                                   (setf highest-codon ac-codon
                                         highest-gc    ac-codon-gc))
                              finally (return (list highest-codon highest-gc)))
          do (unless (or (= ngc 3) ; codon already has 3 Gs and/or Cs
                         ;; codon already has highest number of Gs and/or Cs
                         (= ngc (second highest)))
               (when verbose
                 (format t "~4D: replacing ~A -> ~A (~A)~%"
                         i codon (first highest) ac))
               (setf (elt codons i) (first highest))))
    codons))


;; Same as `replace-with-higher-gc-codon-01` but uses last codon with highest
;; G and/or C count for the aminoacid.  Quite a difference!
(defun replace-with-higher-gc-codon-02 (parsed-side-by-side &key (verbose t))
  (let* ((codons (get-original-codons parsed-side-by-side))
         (tables (make-equivalence-tables))
         (c2a (getf tables :c2a))   ; codon to aminoacid
         (a2c (getf tables :a2c)))  ; aminoacid to codons
    (loop for codon across codons
          for i from 0
          for ac = (gethash codon c2a)      ; get aminoacid for this codon
          for ac-codons = (gethash ac a2c)  ; get all codons for the aminoacid
          for ngc = (n-gc codon)  ; number of Gs and Cs in `codon`
          ;; return last codon with highest number of Gs and Cs
          for highest = (loop with highest-codon = (first ac-codons)
                              with highest-gc = (n-gc highest-codon)
                              for ac-codon in (rest ac-codons)
                              for ac-codon-gc = (n-gc ac-codon)
                              do (when (>= ac-codon-gc highest-gc)
                                   (setf highest-codon ac-codon
                                         highest-gc    ac-codon-gc))
                              finally (return (list highest-codon highest-gc)))
          do (unless (or (= ngc 3) ; codon already has 3 Gs and/or Cs
                         ;; codon already has highest number of Gs and/or Cs
                         (= ngc (second highest)))
               (when verbose
                 (format t "~4D: replacing ~A -> ~A (~A)~%"
                         i codon (first highest) ac))
               (setf (elt codons i) (first highest))))
    codons))


;;; Main

(defun main ()
  (format t "~&~%=== results ===~%")
  (let ((psbs (read-side-by-side-csv))
        (algorithms '(nop
                      berts-algorithm
                      replace-with-higher-gc-codon-01
                      replace-with-higher-gc-codon-02)))
    (loop for algo in algorithms
          for perc = (compare-against-vaccine (funcall algo psbs :verbose nil)
                                              psbs)
          do (format t "~32A: ~,2F%~%" algo perc)))
  (format t "===~%~%You can quit the Lisp environment with `(quit)`.~%"))
