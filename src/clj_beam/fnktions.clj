; Try to outline a Clojure-based version of the Beam Smiles Parsing Library
; https://github.com/johnmay/beam
;

(ns clj-beam.fnktions
  (:require
   [schema.core :as s]
   [plumbing.core :as p :include-macros true]
   [plumbing.fnk.pfnk :as pfnk :include-macros true]
   [plumbing.graph :as graph :include-macros true]
   [plumbing.map :as map])
  ;(:import [uk.ac.ebi.beam Atom Graph Parser CharBuffer])
  )

;;
;;
;; Schemas for Atoms, Bonds , Molecules
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Use the same Graph Model as Beam for holding information about the Chemical Structure
;(defrecord cAtom  isotope element aromatic charge hydrogens atomClass subset)
;(defrecord cGraph atoms edges topologies order size)
;(defrecord cParser stack graph ringbonds arrangement
;                    configurations start openRings strict)


(def Atom
  "An Atom should have an element and optional XYX info"
  {:element [(s/enum :Bromine :Boron :Chlorine :Carbon :Nitrogen
                    :Oxygen  :Phosphorus :Sulfur :Fluorine :Iodine
                    :DEUTERIUM :TRITIUM )]

   (s/optional-key :isotope) long
   (s/optional-key :charge) long
   (s/optional-key :hydrogens) int
   (s/optional-key :x) long
   (s/optional-key :y) long
   (s/optional-key :z) long })

(def Bond
  "A Bond should have an order"
  {:order [(s/enum  :Single :Double :Triple :Quadruple )]
   :aromatic [(s/enum :Yes :No)]
   :type  [(s/enum :Up :Down :Dot)]
   :atoms [Atom]})

(def Molecule
  "A Molecule is a Collection of Atoms and Bonds"
  {(s/optional-key :name) [s/Str]
   :atoms [Atom]
   :bonds  [Bond]})

(def Graph
  "A Graph is a Collection of Atoms and Bonds with some other Information"
  {(s/optional-key :name) [s/Str]
   :atoms  [Atom]
   :bonds  [Bond]
   :topologies  [Bond]
   :order  [Bond]
   :size   long })


;; Prismatic Graph-based Functions for Creating Molecules from Smiles Strings
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;pencil out what these would look like
(def readsmiles
  "A graph specifying the parsing of a smiles string"
  {
   :sqs  (p/fnk [smi] (string->keys smi))
   :atms (p/fnk [sqs] (into [] (map symbol->Atom sqs)))
   :bnds (p/fnk [sqs] (str "needs to be done"))
   :g    (p/fnk [atms] {:atoms atms } )})


(defn string->keys [smi]
  "Take a Smiles String and get a sequence of parsed smiles units"
  (let [conv (fn [x] (cond
                       ( = (apply str x) "Br") :Br
                       ( = (apply str x) "Cl") :Cl
                       :else (keyword (str (first x)))))]
        (->> smi
             (partition 2 1 ); partition into characters
             (filter #(not= (first %) \r)) ;BR abd Cl are special cases
             (filter #(not= (first %) \l))
             (map conv))))


;;
;;
;; Still Need to Figure out the Ring issue

(defn symbol->Atom
  "Map Symbol to Atom Creation"
  [sym]
  (cond
    (in? [:Br :B :Cl :C :N :O :P :S :F :I] sym )
         {:element (sym elementmap) :aromatic :No}
    (in? [:c :n :o :p :s ] sym )
         {:element (sym elementmap) :aromatic :Yes}
    (in? [:H :D :T] sym )
         {:element (sym elementmap) :aromatic :No}
    (in? [:0 :1 :2 :3 :4  :5 :6 :7 :8 :9 ]
         ring)
    :else "There is a prbolem" ))


(def elementmap {
  :* :Unknown    :Br :Bromine   :B  :Boron
  :Cl :Chlorine  :C :Carbon     :N :Nitrogen
  :O :Oxygen     :P :Phosphorus :S :Sulfur
  :F :Fluorine   :I :Iodine     :b :Boron
  :c :Carbon     :n :Nitrogen   :o :Oxygen
  :p :Phosphorus :s :Sulfur     :H :Hydrogen
  :D :DEUTERIUM  :T :TRITIUM })



(def smileseager (graph/compile readsmiles))
(:g  (smileseager {:smi "CCCC"}))


;(def simple "CCCCCC")
;(def smi "C[C]CBr(O)1ccc1CCONC2Br")
;(def smiles-keys (string->keys simple))
;smiles-keys
;(map symbol->Atom smiles-keys)


;(def addatom)
;(def ring)
;(def addaromatic)
;(def elementmap)



;; (def actionmap
;;   "Map Letters to Actions"
;;    {  :Br addatom  :B addatom
;;       :Cl addatom  :C addatom
;;       :N  addatom  :O addatom
;;       :P  addatom  :S addatom
;;       :F  addatom  :I addatom
;;       :b  addaromatic
;;       :c  addaromatic
;;       :n  addaromatic
;;       :o  addaromatic
;;       :p  addaromatic
;;       :s  addaromatic
;;       :H  addatom
;;       :D  addatom
;;       :T  addatom
;;       :0  ring     :1  ring
;;       :2  ring     :3  ring
;;       :4  ring     :5  ring
;;       :6  ring     :7  ring
;;       :8  ring     :9  ring })

(def assemblesmiles
  "A graph specifying how to assemble the process"
  {:smi (p/fnk [smi] (map charmpa ))})

(def assignfunctions
  "A graph specifying which smiles character goes with which function"
  {:smi (p/fnk [smi] (map charmpa ))})


(defn in?
  "true if seq contains elm"
  [seq elm]
  (some #(= elm %) seq))

(defn notin?
  "tru if seq does not contian element"
  [seq elm]
  (not (in? seq elm)))







;(defrecord cParser stack graph ringbonds arrangement
;                    configurations start openRings strict)


;; Addatom
;; addAtom(AtomImpl.AliphaticSubset.Sulfur);
;; addAtom(AtomImpl.AromaticSubset.Boron);
;; g.markDelocalised();
;; case '[':
;;   addAtom(readBracketAtom(buffer));
;; case '9'
;;   ring(c - '0', buffer);
;; case '%':
;;   int num = buffer.getNumber(2);
;;   if (num < 0)
;;       throw new InvalidSmilesException("a number (<digit>+) must follow '%':", buffer);
;;   ring(num, buffer);

;; ;then bonds:
;; case '\\':
;;   bond = Bond.DOWN;

;; // branching
;;     case '(':
;;         if (stack.empty())
;;             throw new InvalidSmilesException("cannot open branch - there were no previous atoms:",
;;                                              buffer);
;;         stack.push(stack.peek());
;;         break;
;;     case ')':
;;         if (stack.size() < 2)
;;             throw new InvalidSmilesException("closing of an unopened branch:",
;;                                              buffer);
;;         stack.pop();
;;         break;

;;     // termination
;;     case '\t':
;;     case ' ':
;;     case '\n':
;;     case '\r':
;;         return;

(defn addatom [at graph]
  "Add an atom to a graph"
  (assoc at g)
  )



;; private void addAtom(Atom a) {
;;         int v = g.addAtom(a);
;;         if (!stack.empty()) {
;;             int u = stack.pop();
;;             if (bond != Bond.DOT)
;;                 g.addEdge(new Edge(u, v, bond));
;;             else
;;                 start.add(v); // start of a new run
;;             if (arrangement.containsKey(u))
;;                 arrangement.get(u).add(v);

;;         }
;;         stack.push(v);
;;         bond = Bond.IMPLICIT;

;;         // configurations used to create topologies after parsing
;;         if (configuration != Configuration.UNKNOWN) {
;;             configurations.put(v, configuration);
;;             configuration = Configuration.UNKNOWN;
;;         }
;;     }

;;     /**


(def funcitonlist {
  :addatom {:members []
            :func ""})

(def charmap {  :* :Unknown
                ;note that disambiguation of the list is needed
                :Br :Bromine
                :B  :Boron
                :Cl :Chlorine
                :C :Carbon
                :N :Nitrogen
                :O :Oxygen
                :P :Phosphorus
                :S :Sulfur
                :F :Fluorine
                :I :Iodine
                ;note that aromatinc is not taken care of
                :b :Boron
                :c :Carbon
                :n :Nitrogen
                :o :Oxygen
                :p :Phosphorus
                :s :Sulfur
                :D :DEUTERIUM
                :T :TRITIUM
                ;note buffer needs dismiguation
                ;:[ :BUFFER
                :0 :RING
                :1 :RING
                :2 :RING
                :3 :RING
                :4 :RING
                :5 :RING
                :6 :RING
                :7 :RING
                :8 :RING
                :9 :RING
                ;:% :temp
                   ; BUFFER.getnumber(2)
                    ;case "%":
                    ;int num = buffer.getNumber(2);
                    ;if (num < 0)
                    ;    throw new InvalidSmilesException("number (<digit>+) must follow "%"", buffer);
                    ;ring(num);
                    ;break;

                ;// bond/dot
                :- :Bond.SINGLE
                := :Bond.DOUBLE
                :# :Bond.TRIPLE
                ;:$ :Bond.QUADRUPLE
                :: :Bond.AROMATIC
                :/ :Bond.UP
                ;:\ :Bond.DOWN
                :. :Bond.DOT
                ;// branching
                ;:( :xx
                    ;if (stack.empty()) throw new InvalidSmilesException("cannot open branch - no previous atom",
                    ;                                     buffer);
                    ;stack.push(stack.peek());
                    ;break;
                ;:) :cvcv
                    ;if (stack.size() < 2)
                    ;    throw new InvalidSmilesException("Closing of an unopened branch",

                ;// termination
                ;:\t :end
                ;:\n :end
                ;:\r :end }}}




;; (for [meth (.getMethods Parser)
;;       :let [name (.getName meth)]]
;;   name)

;; (println .getMethods CharBuffer)
