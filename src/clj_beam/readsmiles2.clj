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
                       :else (keyword (str (first x)))))

        part   (partition 2 1 '(:Padding) smi)   ;partition into characters and use :Padding
        filt1  (filter #(not= (first %) \r) part) ;BR abd Cl are special cases
        filt2  (filter #(not= (first %) \r) filt1)
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
  (let [mapindexed (into [] (map-indexed vector atomlist))
        dblocs     (filter #(= := (first %)) mapindexed)
        tblocs     (filter #(= :# (first %)) mapindexed)
        fb         (filter #(= :( (first %)) mapindexed)
        rb         (filter #(= :) (first %)) mapindexed)
        getsinglebonds (fn [a b] (if (and (map? (second a))
                                          (map? (second b)))
                                     {:order :Single
                                      :aromatic :No
                                      :type  :Dot
                                      :atoms [ (first a) (first b)] }
                                     []
                                   ))

        ;; make the general for all of the bond orders
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        getdoublebonds (fn [x] if (and (atombefore x  mapindexed)
                                       (atomafter  x  mapindexed))
                               (bond :double)
                               [])
        gettriplebonds (fn [x] if (and (atombefore x  mapindexed)
                                       (atomafter  x  mapindexed))
                               (bond :triple)
                               [])

        getnestedbonds  (fn [x] if (and (atombefore x  mapindexed)
                                       (atomafter  x  mapindexed))
                               (bond :triple)
                               [])

        getringbonds   (fn [x] if (and (atombefore x  mapindexed)
                                       (atomafter  x  mapindexed))
                               (bond :triple)
                               [])

        ;; combine bonds
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        singlebonds (filter map?
                (map #(apply getsinglebonds %) (partition 2 1 mapindexed)))

        doublesbonds (filter map?
                (map #(apply getdoublebonds %) (partition 2 1 mapindexed)))

        triplebonds (filter map?
                (map #(apply gettriplebonds %) (partition 2 1 mapindexed)))

        nestedbonds (filter map?
                (map #(apply getnestedbonds %) (partition 2 1 mapindexed)))

        ringbonds (filter map?
                (map #(apply getringbonds %) (partition 2 1 mapindexed)))


                             ]

    (into [] (singlebonds doublebonds triplebonds) )
    ))

(getbonds smitest)


(defn- atombefore [idx vect]
  "for a given index of a vector contianing an indexed vector
  find the first previous vector containing a map/atom "
  (let [v2 (rseq (subvec vect 0 idx))]
       (first (drop-while #(not map? (second %))))))

(defn- atomafter [idx vect]
  "for a given index of a vector containing an indexed vector
  find the next vector contianing a map/atom "
  (let [v2 (subvec vect 0 idx)]
       (first (drop-while #(not map? (second %))))))

(defn- makebond [bondtype pos vect]
  "create bond"
   (let [ab (atombefore pos vect)
         aa (atomafter  pos vect)]
     ()))

; bonds are all explicit double and triple bonds as well as implicit
;

(defn- singlebonds [atomvect]
  "Parse Sequence and Generate Sequence Bonds. Note: This will NOT currently work on anyhting inside of brackets"
   (let [ part  (partition 2 1 atomvect)
          bonds (filter #(and (map? (second (first  %)))
                              (map? (second (second %)))) part)
          singlebond (fn [bond]
                       {:atom1 (first (first bond)) :atom2 (first (second bond)) :order :Single :aromatic :No})]
          (into [] (map singlebond bonds))))

(defn- rings    [atomvect]
  "Bonds from Rings"
   (let [nums  (filter #(integer? (second  %)) atomvect) ;getrings numbers

         ]
     ;add rings based on numbers check the
         "a"))

(defn- doublebonds  [atomvect]
  "Detect Double Bonds: Note that this method does not see things in brackets"
     (let [b  (filter #(=  := (second  %)) atomvect) ]
       ;check for being left of a parenthesis
       ))

(defn- triplebonds  [atomvect]
  "Detect Triple Bonds: Note that this method doe snot see things in brackets"
     (let [b  (filter #(=  := (second  %)) atomvect) ]
       ;check for being left of a parenthesis
       ))

(defn- bondparenthesis  [atomvect]
  "Detect bonds across parentheses"
     (let []
       ;check all nested and long bond
       ))



(def smiles (string/split  (slurp "data/smiles.txt") #"\n"))
(def smi12 (nth smiles 12))
(map readsmiles smiles)
smiles
smi12
(def smitest (readsmiles smi12))
smitest
(getbonds smitest)
