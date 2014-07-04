(ns clj-beam.readsmiles
  (:gen-class)
  (:require
   [clojure.set :as cset]
   [clojure.string :as string]
   [schema.core :as s]
   [plumbing.core :as p :include-macros true]
   [plumbing.fnk.pfnk :as pfnk :include-macros true]
   [plumbing.graph :as graph :include-macros true]
   [plumbing.map :as map]
;   [clj-beam.schemas :as schemas]
   [clj-beam.regex :only bracketatominfo]
;   [clj-beam.opensmiles :only element_symbols]
   ))


(defn string->Atom [smi]
  "take a smiles string and return a list of atoms or keyword characters

    (string->Atom 'CCc1')
       ( {:element :Carbon, :aromatic :No, :index 0}
         :(
         {:element :Carbon, :aromatic :No, :index 2}
         :)
         {:element :Oxygen, :aromatic :No, :index 4})"
  (let [part   (partition 2 1 '(:Padding) smi)   ;partition into characters and use :Padding
        filt1  (filter #(not= (first %) \r) part) ;BR abd Cl are special cases
        filt2  (filter #(not= (first %) \r) filt1)
        smivec (vec (map conv filt2))
        indx   (range (count smivec))
        in?    (fn [seq elm] (some #(= elm %) seq))
        conv   (fn [x] (cond
                       ( = (apply str x) "Br") :Br
                       ( = (apply str x) "Cl") :Cl
                       :else (keyword (str (first x))))) ]

    (for [ idx indx ]
       (let [ sym (get smivec idx)]
         (cond
          (in? [:Br :B :Cl :C :N :O :P :S :F :I] sym )
               {:element ( sym elementmap) :aromatic :No :index idx}
          (in? [:c :n :o :p :s ] sym )
               {:element (sym elementmap) :aromatic :Yes}
          (in? [:H :D :T] sym )
               {:element (sym elementmap) :aromatic :No}
          :else sym )))))


(defn bracketatominfo [a]
  "a bracketed   atom can may have the following:
      isotope:    a nubmer before the symbol
      symbol:     the atomic symbol
      chirality:  for chiral centers
      protons:    H with an optional digit after
      charge:     multiple +/- with optional digit
      class:      optional field to specify extra information
  return a map with atom information"
  (let [;get basic charge, proton, and isotope information
        element   (let  [ el (re-seq #"(?i)[abcdefgijklmnopqrstuvmwyz]+" a) ]
                    (if (nil? el)  "H" el)) ;pattern excludes "H" for ease. add it back if necessary
        isotope   (re-find #"^[0-9]+" a)
        chiral    (re-seq #"[@]+" a)
        chiralH   (re-seq #"[@]+[H]" a)
        charge    (first (re-seq #"[+-]+" a))
        chargenum (first (map read-string (re-find #"[+-](\d)" a)))
        protons   (re-seq #"[H]+" a)
        protnum   (second ( map read-string (re-find #"[H](\d)" a)))
        atomclass (last (re-seq #":(?i)(\w+)" a))

        ;resolve ambiguities with charge and protons. this could probably be handled
        ;better with better knowldege of regex.
        fin_charge (cond
                     chargenum chargenum
                     charge  ( if (= (first charge) \+)
                                  (count charge)
                                  (- (count charge)))
                     :else 0)

        fin_protons (cond
                       protnum protnum
                       protons 1
                       :else 0)]

    {  :element   (first element)
       :isotope   isotope
       :charge    fin_charge
       :hydrogens fin_protons
       :atomclass atomclass }
    ))


(defn readsmiles [smi]
  "Parse a smiles string and return a list that contains atoms and bond-symbols"
  (let [startswithbracket (if (= \[ (first smi)) true false)
        brackets (re-pattern #"(\[)|(\])") ; bracket pattern
        splitsmi (partition-all 2 (string/split smi brackets) ) ;splits on brackets every other sequence was bracketed
        f  (map first splitsmi)
        s  (remove nil? (map second splitsmi))
        inter (fn [a b] (if (= (count a) (count b))
                            (interleave a b)
                            (cons (first a) (interleave b (rest a)))))]
     (if startswithbracket
       (let [fm (map bracketatominfo f)
             rm (map string->Atom s)]
                (flatten (inter fm rm )))
       (let [fm (map string->Atom f)
             rm (map bracketatominfo s)]
                (flatten (inter fm rm ))))))

(def smiles (string/split  (slurp "data/smiles.txt") #"\n"))
(def smi12 (nth smiles 12))
(map readsmiles2 smiles)