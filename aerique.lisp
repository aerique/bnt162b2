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


;;; Globals

(defparameter *a2c* nil)  ; maps aminoacid to its codon(s)
(defparameter *c2a* nil)  ; maps codon to the aminoacid it creates


;;; Common Functions

(defun mk-pathname (pathname-or-string)
  (if (typep pathname-or-string 'string)
      (merge-pathnames pathname-or-string)
      pathname-or-string))


(defun print-hash-table (hash-table)
  (maphash (lambda (k v)
             (format t "~S: ~S~%" k v))
           hash-table))


(defun random-elt (sequence)
  (elt sequence (random (length sequence))))


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

;; FIXME don't do this here but in Globals at the top
(let ((equivalence-tables (make-equivalence-tables)))
  (setf *a2c* (getf equivalence-tables :a2c)
        *c2a* (getf equivalence-tables :c2a)))


;;; Functions

(defun codons-to-nucleotides (codons)
  (loop for codon across codons
        append (list (elt codon 1) (elt codon 1) (elt codon 2)) into result
        finally (return (coerce result 'string))))


(defun compare-codons-against-vaccine (codons parsed-side-by-side)
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


(defun compare-nucleotides-against-vaccine (codons parsed-side-by-side)
  "CODONS should be a vector of codons where each codon is a three letter
  string: #(\"TAA\" \"ACA\" ...)
  Returns the percentage of nucelotides EQUAL to the ones in the vaccine."
  (loop with vaccine-nucleotides = (codons-to-nucleotides
                                    (get-vaccine-codons parsed-side-by-side))
        with n_different = 0
        for nucleotide across (codons-to-nucleotides codons)
        for i from 0
        do (when (char/= nucleotide (elt vaccine-nucleotides i))
             (incf n_different))
        finally (return (- 100 (* (/ n_different (length vaccine-nucleotides))
                                  100.0)))))


(defun get-virus-codons (parsed-side-by-side)
  "PARSED-SIDE-BY-SIDE is the output of READ-SIDE-BY-SIDE-CSV."
  (loop for row across parsed-side-by-side
        collect (elt row 1) into result
        finally (return (coerce result 'vector))))


(defun get-vaccine-codons (parsed-side-by-side)
  "PARSED-SIDE-BY-SIDE is the output of READ-SIDE-BY-SIDE-CSV."
  (loop for row across parsed-side-by-side
        collect (elt row 2) into result
        finally (return (coerce result 'vector))))


;;; No Operation 'Algorithm

(defun nop (parsed-side-by-side &key verbose)
  (declare (ignore verbose))
  (get-virus-codons parsed-side-by-side))


;;; Bert's Sample Algorithm: ../00-Resources/bnt162b2-git/3rd-gc.go

(defun berts-algorithm (parsed-side-by-side &key (verbose t))
  (let ((codons (get-virus-codons parsed-side-by-side))
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
  (let* ((codons (get-virus-codons parsed-side-by-side))
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
  (let* ((codons (get-virus-codons parsed-side-by-side))
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


(defun replace-cca-with-cct (parsed-side-by-side &key (verbose t))
  (let ((codons (get-virus-codons parsed-side-by-side)))
    (loop for codon across codons
          for i from 0
          do (when (string= codon "CCA")
               (when verbose
                 (format t "~4D: replacing CCA -> CCT~%" i))
               (setf (elt codons i) "CCT")))
    codons))


;;; GA Static Remap Attempt

(defun make-random-codon-mapping ()
  (let ((mapping (make-hash-table :size 64 :test #'equal)))
    (loop for codon in (loop for codon being the hash-keys in *c2a*
                             collect codon)
          for aminoacid = (gethash codon *c2a*)
          for codons = (gethash aminoacid *a2c*)
          do (setf (gethash codon mapping) (random-elt codons)))
    mapping))  ; hash table (faster than an alist (hopefully (pretty sure)))

(defun make-random-codon-mapping-test ()
  (let ((mapping (make-hash-table :size 64 :test #'equal)))
    (loop for codon in (loop for codon being the hash-keys in *c2a*
                             repeat 8
                             collect codon)
          for aminoacid = (gethash codon *c2a*)
          for codons = (gethash aminoacid *a2c*)
          do (setf (gethash codon mapping) (random-elt codons)))
    mapping))  ; hash table (faster than an alist (hopefully (pretty sure)))


(defun copy-codon-mapping (codon-mapping)
  (let ((copy (make-hash-table :size 64 :test #'equal)))
    (maphash (lambda (k v)
               (setf (gethash k copy) v))
             codon-mapping)
    copy))


(defun print-codon-mapping (codon-mapping)
  "CODON-MAPPING is the output of MAKE-RANDOM-CODON-MAPPING (a hash table)."
  (maphash (lambda (k v)
             (format t "[~A->~A] " k v))
           codon-mapping)
      (format t "~%"))


(defun do-remap (virus-codons codon-mapping)
  (loop for codon across virus-codons
        collect (gethash codon codon-mapping) into result
        finally (return (coerce result 'vector))))


(defun ga-static-remap (parsed-side-by-side &key (codon-mapping nil)
                                                 (generations 10000)
                                                 (verbose t))
  (let ((mapping (if codon-mapping
                     (copy-codon-mapping codon-mapping)
                     (make-random-codon-mapping)))
        (virus-codons (get-virus-codons parsed-side-by-side)))
    ;(when verbose
    ;  (format t "mapping: ")
    ;  (print-codon-mapping mapping))
    (loop with best-mapping          = nil
          with best-codon-match      = 0
          with best-nucleotide-match = 0
          repeat generations
          for i from 0
          for new-codons           = (do-remap virus-codons mapping)
          for new-codon-match      = (compare-codons-against-vaccine
                                      new-codons parsed-side-by-side)
          for new-nucleotide-match = (compare-nucleotides-against-vaccine
                                      new-codons parsed-side-by-side)
          do ;; Note: If the matches are exactly the same the new codons are
             ;;       skipped, this might ignore codons that provide a better
             ;;       folding match!
             ;;
             ;; A codon match is preferred over a nucleotide match.  This is
             ;; an arbitrary choice.
             (when (or (> new-codon-match best-codon-match)
                       (and (= new-codon-match      best-codon-match)
                            (> new-nucleotide-match best-nucleotide-match)))
               (setf best-mapping          (copy-codon-mapping mapping)
                     best-codon-match      new-codon-match
                     best-nucleotide-match new-nucleotide-match)
               (when verbose
                 (format t "  - [~8D] improved generation: codon-match=~2,2F% ~
                            nucleotide-match=~2,2F%~%"
                         i best-codon-match best-nucleotide-match)
                 (force-output)))
             ;; Change a random mapping.
             ;; FIXME call `mutate` now that we have it
             (let* ((keys (loop for k being the hash-keys in mapping
                                collect k))
                    (random-key (random-elt keys))
                    ;; We use one of the mapping codons (`random-key`) to get
                    ;; its aminoacid and use that aminoacid the get the
                    ;; equivalent codons and pick a random one of those.
                    (new-codon (random-elt (gethash (gethash random-key *c2a*)
                                                    *a2c*))))
               ;(when verbose
               ;  (format t "  - [~8D] replacing ~A -> ~A~%"
               ;          i (gethash random-key mapping) new-codon)
               ;  (force-output))
               (setf (gethash random-key mapping) new-codon))
          finally (return best-mapping))))


;; Experience so far is that in the initial 10000 to 100000 generations the
;; most optimizations are found and then the GA run gets stuck in a local
;; optimum.  So we do multiple shorter runs and see how that goes.  Then we
;; can spend more time on the best result.
(defun multiple-ga-runs (parsed-side-by-side &key (generations 1000)
                                                  (runs 100) (verbose t))
  (let ((virus-codons (get-virus-codons parsed-side-by-side)))
    (when verbose
      (format t "=== doing ~D GA runs ===~%" runs)
      (force-output))
    (loop with best-mapping          = nil
          with best-codon-match      = 0
          with best-nucleotide-match = 0
          repeat runs
          for i from 0
          for new-mapping          = (ga-static-remap parsed-side-by-side
                                                      :generations generations
                                                      :verbose nil)
          for new-codons           = (do-remap virus-codons new-mapping)
          for new-codon-match      = (compare-codons-against-vaccine
                                      new-codons parsed-side-by-side)
          for new-nucleotide-match = (compare-nucleotides-against-vaccine
                                      new-codons parsed-side-by-side)
          do ;(when verbose
             ;  (format t "--- run ~D ---~%" i))
             ;; See `ga-static-remap`.
             (when (or (> new-codon-match best-codon-match)
                       (and (= new-codon-match      best-codon-match)
                            (> new-nucleotide-match best-nucleotide-match)))
               (setf best-mapping          (copy-codon-mapping new-mapping)
                     best-codon-match      new-codon-match
                     best-nucleotide-match new-nucleotide-match)
               (when verbose
                 (format t "o [~6D] improved mapping: codon-match=~2,2F% ~
                            nucleotide-match=~2,2F%~%"
                         i best-codon-match best-nucleotide-match)
                 (force-output)))
          finally (return best-mapping))))


(defun make-population (&key (size 100) (verbose t))
  (declare (ignore verbose))
  (loop for i from 0 below size
        collect (make-random-codon-mapping) into result
        finally (return (coerce result 'vector))))


(defun determine-population-fitness (population parsed-side-by-side)
  (loop with psbs = parsed-side-by-side
        with virus-codons = (get-virus-codons parsed-side-by-side)
        for individual across population
        for codons = (do-remap virus-codons individual)
        for codon-match = (compare-codons-against-vaccine codons psbs)
        for nucleotide-match = (compare-nucleotides-against-vaccine codons psbs)
        for fitness = (+ (* 2 codon-match) nucleotide-match)
        collect (cons fitness individual) into result
        finally (return (sort (coerce result 'vector) #'> :key #'car))))


(defun tournament-selection (fitness-population &key (rounds 4))
  "FITNESS-POPULATION is the output of DETERMINE-POPULATION-FITNESS."
  (loop with best = (random-elt fitness-population)
        repeat rounds
        for contestant = (random-elt fitness-population)
        do (when (> (car contestant) (car best))
             (setf best contestant))
        finally (return (cdr best))))


(defun mapping-to-vector (mapping)
  (loop for k being the hash-keys in mapping
        ;collect (vector k (gethash k mapping)) into result
        ;; alternative visualization
        ;collect (list k '-> (gethash k mapping)) into result
        collect (cons k (gethash k mapping)) into result
        finally (return (coerce result 'vector))))


(defun mate (mapping-a mapping-b)
  (let* ((mv-a (mapping-to-vector mapping-a))
         (mv-b (mapping-to-vector mapping-b))
         (pivot (random (length mv-a)))  ; same length as `mv-b`
         (mv-c (concatenate 'vector (subseq mv-a 0 pivot)
                                    (subseq mv-b pivot)))
         (mv-d (concatenate 'vector (subseq mv-b 0 pivot)
                                    (subseq mv-a pivot)))
         (mapping-c (make-hash-table :size 64 :test #'equal))
         (mapping-d (make-hash-table :size 64 :test #'equal)))
    (loop for cons across mv-c
          for from = (car cons)
          for to   = (cdr cons)
          do (setf (gethash from mapping-c) to))
    (loop for cons across mv-d
          for from = (car cons)
          for to   = (cdr cons)
          do (setf (gethash from mapping-d) to))
    (list mapping-c mapping-d)))


(defun mutate (mapping)
  ;; Change a random mapping.
  (let* ((new-mapping (copy-codon-mapping mapping))
         (keys (loop for k being the hash-keys in mapping collect k))
         (random-key (random-elt keys))
         ;; We use one of the mapping codons (`random-key`) to get
         ;; its aminoacid and use that aminoacid the get the
         ;; equivalent codons and pick a random one of those codons.
         (new-codon (random-elt (gethash (gethash random-key *c2a*) *a2c*))))
    ;(format t "Changing ~A->~A to ~A->~A.~%"
    ;        random-key (gethash random-key mapping) random-key new-codon)
    (setf (gethash random-key new-mapping) new-codon)
    new-mapping))


(defun evolve (population parse-side-by-side)
  (loop with size = (length population)
        with f-p = (determine-population-fitness population parse-side-by-side)
        with new-population = ;; elitisism: fill part of new population with
                              ;; the best member of the previous population
                              (loop for i across (subseq f-p 0 (floor size 10))
                                    collect (cdr i))
        while (< (length new-population) (length population))
        do (if (<= (random 1.0) 0.8)
               ;; crossover
               (loop for child in (mate (tournament-selection f-p)
                                        (tournament-selection f-p))
                     do (if (<= (random 1.0) 0.03)
                            (push (mutate child) new-population)
                            (push child new-population)))
               ;; no crossover, no mating
               (if (<= (random 1.0) 0.03)
                   (push (mutate (random-elt population)) new-population)
                   (push (random-elt population) new-population)))
        finally (return (coerce new-population 'vector))))


(defun do-runs (parsed-side-by-side
                &key (max-generations 10) (population-size 16))
  (let* ((p  (make-population :size population-size))
         (pf (determine-population-fitness p parsed-side-by-side)))
    (format t "[~6D] best: ~S~%" 0 (elt pf 0))
    (force-output)
    (loop repeat max-generations
          for i from 1
          do (setf p  (evolve p parsed-side-by-side)
                   pf (determine-population-fitness p parsed-side-by-side))
             (format t "[~6D] best: ~S~%" i (elt pf 0))
             (force-output))
    p))


;;; Main

(defun main ()
  (format t "~&~%=== results ===~%~%")
  (let ((psbs (read-side-by-side-csv))
        (algorithms '(nop
                      berts-algorithm
                      replace-with-higher-gc-codon-01
                      replace-with-higher-gc-codon-02
                      replace-cca-with-cct)))
    (format t " name                             | codon match | nucleotide match~%~
               ----------------------------------+-------------+------------------~%")
    (loop for algo in algorithms
          for codons = (funcall algo psbs :verbose nil)
          for codon_perc = (compare-codons-against-vaccine codons psbs)
          for nucl_perc = (compare-nucleotides-against-vaccine codons psbs)
          do (format t " ~32A |      ~2,2F% | ~2,2F%~%"
                     algo codon_perc nucl_perc)))
  (format t "~%===~%~%You can quit the Lisp environment with `(quit)`.~%"))
