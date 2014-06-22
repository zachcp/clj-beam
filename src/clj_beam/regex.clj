(ns clj-beam.reagex
  (:require
   [clojure.set :as cset]
   [clojure.string :as string]
   [schema.core :as s]
   [clj-beam.schemas :as schemas]))


; example Strings to pick out
; OH3+ 2H 235U C@@H C@H I- Na+ Br- H+ Fe+2 OH- Fe++ OH3+ NH4+

(def examples ["OH3+"  "2H" "235U" "C@@H" "C@H" "I-" "Na+" "Br-" "H+" "Fe+2" "OH-" "Fe++" "OH3+" "NH4+"])

(defn exampleatom [a]
  (let [;get basic charge, proton, and isotope information
        charge    (re-seq #"[+-]+" a)       ;get charge
        chargenum (re-seq #"[+-]+[0-9]" a)  ;get charge if specified by number
        protons   (re-seq #"[H]+" a)
        protnum   (re-seq #"[H][0-9]" a)
        isotope   (re-seq #"(?i)[0-9]+[abcdefghijklmnopqrstuvmwyz]" a)
        ;check chirality
        chiral   (re-seq #"[@]+" a)
        ]
    [a charge chargenum  protons protnum  isotope chiral]))

(map exampleatom examples)

