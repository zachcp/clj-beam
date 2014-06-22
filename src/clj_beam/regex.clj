(ns clj-beam.reagex
  (:require
   [clojure.set :as cset]
   [clojure.string :as string]
   [schema.core :as s]
   [clj-beam.schemas :as schemas]))


; example Strings to pick out
; OH3+ 2H 235U C@@H C@H I- Na+ Br- H+ Fe+2 OH- Fe++ OH3+ NH4+

(def examples ["OH3+"  "35Cl:classinfo" "2H" "235U" "C@@H" "C@H" "I-" "Na+" "Br-" "H+" "Fe+2" "OH-" "Fe++" "OH3+" "NH4+"])

(defn bracketatominfo [a]
  "a bracketed   atom can may have the followinwg:
     isotope:    a nubmer before the symbol
     symbol:     the atomic symbol
     chirality:  for chiral centers
     protons:    H with an optional digit after
     charge:     multiple +/- with optional digit
     class:      optional field to specify extra information"
  (let [;get basic charge, proton, and isotope information
        element   (re-seq #"(?i)[abcdefghijklmnopqrstuvmwyz]+" a)
        isotope   (re-seq #"(?i)[0-9]+[abcdefghijklmnopqrstuvmwyz]" a)
        chiral    (re-seq #"[@]+" a)
        chiralH   (re-seq #"[@]+[H]" a)
        charge    (re-seq #"[+-]+" a)
        chargenum (re-seq #"[+-]+[0-9]" a)
        protons   (re-seq #"[H]+" a)
        protnum   (re-seq #"[H][0-9]" a)
        atomclass (re-seq #"[:](?i)(\w+)" a)
        tempatom  {}]
    [a element isotope chiral chiralH charge chargenum protons protnum atomclass]))


(map bracketatominfo examples)

