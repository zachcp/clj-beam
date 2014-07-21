(ns clj-beam.readsmiles2
  (:gen-class)
  (:require
   [clojure.set :as cset]
   [clojure.string :as string]
   [schema.core :as s]))


(def elementmap {
  :* :Unknown    :Br :Bromine   :B  :Boron
  :Cl :Chlorine  :C :Carbon     :N :Nitrogen
  :O :Oxygen     :P :Phosphorus :S :Sulfur
  :F :Fluorine   :I :Iodine     :b :Boron
  :c :Carbon     :n :Nitrogen   :o :Oxygen
  :p :Phosphorus :s :Sulfur     :H :Hydrogen
  :D :DEUTERIUM  :T :TRITIUM })

(defn- string->Atom [smi]
  "take a smiles string and return a list of atoms or keyword characters
    (string->Atom 'CCc1')
       ( {:element :Carbon, :aromatic :No, :index 0}
         :(
         {:element :Carbon, :aromatic :No, :index 2}
         :)
         {:element :Oxygen, :aromatic :No, :index 4})"
  (let [in?    (fn [seq elm] (some #(= elm %) seq))
        conv   (fn [x] (cond
                       ( = (apply str x) "Br") :Br
                       ( = (apply str x) "Cl") :Cl
                       ( re-find #"\d" (str (first x))) (Integer/parseInt (str (first x )))
                       :else (keyword (str (first x)))))

        part   (partition 2 1 '(:Padding) smi)   ;partition into characters and use :Padding
        filt1  (filter #(not= (first %) \r) part) ;BR abd Cl are special cases
        filt2  (filter #(not= (first %) \l) filt1)
        smivec (vec (map conv filt2))
        indx   (range (count smivec))]

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

(defn- bracketatominfo [a]
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

    (do
      (println chiral chiralH))
    {  :element   (first element)
       :isotope   isotope
       :charge    fin_charge
       :hydrogens fin_protons
       :atomclass atomclass
       :chiral    (if chiralH (first chiralH) (first chiral))}
    ))


(defn readsmiles [smi]
  "Parse a smiles string and return a list that contains atoms and bond-symbols"
  (let [startswithbracket (if (= \[ (first smi)) true false)
        brackets (re-pattern #"(\[)|(\])") ; bracket pattern
        splitsmi (partition-all 2 (string/split smi brackets) ) ;splits on brackets every other sequence was bracketed
        f  (map first splitsmi)
        s  (remove nil? (map second splitsmi))
        inter (fn [a b]
                (if (= (count a) (count b))
                    (interleave a b)
                    (cons (first a) (interleave b (rest a)))))]
     (if startswithbracket
       (let [fm (map bracketatominfo f)
             rm (map string->Atom s)]
                (into [] (flatten (inter fm rm ))))
       (let [fm (map string->Atom f)
             rm (map bracketatominfo s)]
                (into [] (flatten (inter fm rm )))))))

(defn- getbonds [atomlist]
  "Bonds From Atom List"
  (let [singlebonds (getsinglebonds atomlist)
        doublebonds (getdoublebonds atomlist)
        triplebonds (gettriplebonds atomlist)
        ringbonds   (getbrokenringbonds atomlist)
        ringbonds2  (getadjacentringbonds atomlist)]
    (into [] (concat singlebonds doublebonds triplebonds ringbonds ringbonds2))))

(defn- getindices [x coll]
  "return the indeces of the collection whose value is x"
  (map first
     (filter #(= (second %) x)
             (map-indexed vector coll))))

(defn- getsinglebonds [atomlist]
  "returns single bonds from adjacent atoms"
  (let [bonds (for [i (range (count atomlist))]
                (cond (and (map? (nth atomlist i))
                           (map? (nth atomlist (+ i 1))))
                       {:order :Single
                        :aromatic :No
                        :type  :Dot
                        :atoms [ i (+ i 1)]}))]
    (filter #(not (nil? %)) bonds)))

(defn- getdoublebonds [atomlist]
  "returns double bonds from adjacent atoms
   must take into account:
        direct:   C=C
        indirect: C(=O)
        nested:   C(C)(=O)"
  (let [dblocs   (getindices := atomlist)
        checkmap (fn [x] (map? (nth x atomlist)))
        direct?  (fn [x] (if (and (checkmap (inc x))
                                 (checkmap (dec x)))
                             true false ))
        indirect? (fn [x] (if (and (checkmap (inc x))
                                  (checkmap (dec (dec x))))
                             true false ))]
    (for [d dblocs]
      (cond
         (direct? d)   {:order :Double :atoms [(dec d) (inc d)] }
         (indirect? d) {:order :Double :atoms [ (getnestedatomindex (dec d) coll 0 0) (inc d)]}))))

(defn- gettriplebonds [atomlist]
  "returns double bonds from adjacent atoms
   must take into account:
        direct:   C#C
        indirect: C(#O)
        nested:   C(C)(#O)"
  (let [tblocs   (getindices :# atomlist)
        checkmap (fn [x] (map? (nth x atomlist)))
        direct?  (fn [x] (if (and (checkmap (inc x))
                                 (checkmap (dec x)))
                             true false ))
        indirect? (fn [x] (if (and (checkmap (inc x))
                                  (checkmap (dec (dec x))))
                             true false ))]
    (for [d tblocs]
      (cond
         (direct? d)   {:order :Double :atoms [(dec d) (inc d)] }
         (indirect? d) {:order :Double :atoms [ (getnestedatomindex (dec d) atomlist 0 0) (inc d)]}))))

(defn- getnestedatomindex [idx coll lb rb]
  "return the atom location for an implicit, nested bond
   examples:
     CCN(CC)CC
     CC(=O)O
     NC(C)C(=O)O
     CC(C)(C(=O)O)Oc1ccc(cc1)Cl
     COc1cc(c(c2c1OCO2)OC)CC=C
     Cc1ccccc1NC(=O)C(C)N2CCCC2
     CC(=O)Oc1ccccc1C(=O)[O-]
     c1cc(ccc1C(CC(=O)O)CN)Cl
  x should initially be the position immedietly to the right of the closed right bracket"
;;   (do
;;     (println idx lb rb)
;;     (println (= lb rb))
;;     (println (= :( (nth coll idx) ))
;;     (println (= :) (nth coll idx) )))
  (cond
     (< idx 0)
       "Somethingwrong"
     (and (= lb rb) (> rb 0))
       idx
     (= :( (nth coll idx))
       (getnestedatomindex (dec idx) coll (inc lb) rb)
     (= :) (nth coll idx))
       (getnestedatomindex (dec idx) coll lb (inc rb))
     :else (getnestedatomindex (dec idx) coll lb rb)))

(defn- getbrokenringbonds [atomlist]
  "Return Bonds from Broken rings and Single Bonds
   For now just return the atom immedietely to the left of the integer.
  ToDO: add use cases if there are multiple integers associated with the same atom"
  (let [mapindexed (into [] (map-indexed vector atomlist))
        ringlocs   (filter #(integer? (second %)) mapindexed)
        matches    (for [x ringlocs y ringlocs
                        :when (= (second x) (second y) )
                        :while (not= (first x) (first y))]
                        [x y]) ]
      (for [[a b] matches]
            (let [mi (min (first a) (first b))
                  ma (max (first a) (first b))]
               {:order :Single
                :aromatic :No
                :type  :Dot
                :atoms [ (dec mi) (dec ma)]}
              ))))


(defn- getadjacentringbonds [atomlist]
  "ring integers break up the detection of single rings find ring locations and add the bonds"
   (let [mapindexed (into [] (map-indexed vector atomlist))
         ringlocs   (filter #(integer? (second %)) mapindexed)]
;     ringlocs))
     (for [ [i n] ringlocs]
        (cond
            (>= i (- (count atomlist) 1))
               nil
            (and (map? (nth atomlist (inc i))) (map? (nth atomlist (dec i))))
               {:order :Single
                :aromatic :No
                :type  :Dot
                :atoms [ (dec i) (inc i )]}
             :else "Problem!"  ))))


(def smiles (string/split  (slurp "data/smiles.txt") #"\n"))
(def smi12 (nth smiles 12))
(def smiA (nth smiles 6))
smiA
smiles

(def a (map readsmiles smiles))

(getbonds (nth (filter vector? a) 7))
a
smiles
smi12
smitest
(getbonds smitest)

(def a (map read-string "c1cccc1"))
(integer? (second (map read-string (map str (seq "C1ccc1")))))
