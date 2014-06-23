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

(def Arangement
  "Arrangements are used to specify the atoms that each atom is connected to"
  { [s/Int] [s/Int] })



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




;; ;; Near Literal Translation of the Parser
;; (defrecord parser [stack graph ringbonds arrangement configurations start bond openrings strict])
;; (defrecord graph [atoms edges topologies order size delocalized])

;; (defn baseparser []
;;   "Starting Point for Parsing"
;;    (parser.
;;     (into [] ( range 10)) ;stack
;;     (graph. [] ;atoms
;;             [] ;edges
;;             [] ;topologies
;;             0  ;order;size
;;             1
;;             :Yes) ;graph

;;     ))

(def baseparser
  {:stack '(1 2 3 4 5 6 7 8 9 10)
   :graph basegraph
   :ringbonds [] ;list of rings
   :arrangement {}
   :configurations "none"
   :start #()
   :bond :IMPLICIT
   :openRings 0
   :strict nil})

(def basegraph
  { :atoms  []
    :edges  []
    :topologies  {}
    :order  0
    :size   0
    :delocalized :True})

(pop  '(1 2 3 4 5 6 7 8 9 10))
(peek '(1 2 3 4 5 6 7 8 9 10))
(pop (range))


(defn addatom [symb parser idx atomtype]
  ;parser.Java line 265
  ;grapho.Java line 108

  ;add atom
  (assoc-in parser [:graph :atoms idx ] {:element (symb elementmap) :aromatic atomtype} )

  ;get bond information and updae bond and association information
  (let [ v ( - (:order(:graph parser)) 1) ;order
         l (count (:stack parser))
         u (peek  (:stack parser))]

    (println v l u)
    ;add bond
    (when (not= (:bond parser) :Bond.DOT) ;add a bond
          (assoc-in parser [:graph :edge idx] [u v (:bond parser)] ))
    ;update arrangements
    (when (in? (keys (:arrangement parser)) u)
          (let [c (count (get (:arrangement parser) u))]
             (assoc-in parser [:arrangement u (inc c) ] v )
            ))
    ;Update the stack
    (let [stack (:stack parser)]
      (assoc parser :stack (cons v stack)))
    ;reset the bond type
    (assoc parser :bond :bond.IMPLICIT))

  parser)

 (def a  (addatom :Br baseparser 0 :aromatic))
 (def b  (addatom :Br a 0 :aromatic))

 b

a
(assoc-in a
(cons 1 '(2 3 4))
   private LocalArrangement createArrangement(int u) {
        LocalArrangement la = arrangement.get(u);
        if (la == null) {
            la = new LocalArrangement();
            for (Edge e : g.edges(stack.peek()))
                la.add(e.other(u));
            arrangement.put(u, la);
        }
        return la;
    }


      private void addAtom(Atom a) {
        int v = g.addAtom(a);

if (!stack.empty()) {
            int u = stack.pop();
            if (bond != Bond.DOT)
                g.addEdge(new Edge(u, v, bond));
            else
                start.add(v); // start of a new run
            if (arrangement.containsKey(u))
                arrangement.get(u).add(v);

        }
        stack.push(v);
        bond = Bond.IMPLICIT;

        // configurations used to create topologies after parsing
        if (configuration != Configuration.UNKNOWN) {
            configurations.put(v, configuration);
            configuration = Configuration.UNKNOWN;
        }
    }



(defn parsesmiles [parser symb]
  ;try a literal translation ofhte parsing function
  (let [idx (count (:atoms (:graph parser)))]
    (cond
      (in? [:Br :B :Cl :C :N :O :P :S :F :I :H :D :T] symb )
          (;add an atom to the atomlist
           (assoc-in parser [:graph :atoms (inc idx)] {:element (symb elementmap) :aromatic :No})
           (addatom parser idx atomtype)
           ;add an edge
           (let [v (:order(:graph parser))
                 u (pop (:stack parser))
                 b-idx (count (:bonds (:graoh parser)))]
              (assoc-in parser [:graph :bonds (inc b-idx)] {:element (symb elementmap) :aromatic :No})
              (if (u in (keys (:arrangement parser)))
               ;line 274 of the parser program
                assoc-in parser [:graph :bonds (inc b-idx)] {:element (symb elementmap) :aromatic :No})))
      (in? [:c :n :o :p :s ] symb )
          (assoc-in parser [:graph :atoms (inc idx)] {:element (symb elementmap) :aromatic :Yes}  )
      (in? [:0 :1 :2 :3 :4  :5 :6 :7 :8 :9 ]
         ring)

     (= :* symb)
     (= :B symb)
          (assoc-in parser [:graph :atoms (inc idx)] {:element :Boron}  )


    )))
                case 'B':
                    if (buffer.getIf('r'))
                        addAtom(AtomImpl.AliphaticSubset.Bromine);
                    else
                        addAtom(AtomImpl.AliphaticSubset.Boron);
                    break;
                case 'C':
                    if (buffer.getIf('l'))
                        addAtom(AtomImpl.AliphaticSubset.Chlorine);
                    else
                        addAtom(AtomImpl.AliphaticSubset.Carbon);
                    break;
                case 'N':
                    addAtom(AtomImpl.AliphaticSubset.Nitrogen);
                    break;
                case 'O':
                    addAtom(AtomImpl.AliphaticSubset.Oxygen);
                    break;
                case 'P':
                    addAtom(AtomImpl.AliphaticSubset.Phosphorus);
                    break;
                case 'S':
                    addAtom(AtomImpl.AliphaticSubset.Sulfur);
                    break;
                case 'F':
                    addAtom(AtomImpl.AliphaticSubset.Fluorine);
                    break;
                case 'I':
                    addAtom(AtomImpl.AliphaticSubset.Iodine);
                    break;

                // aromatic subset
                case 'b':
                    addAtom(AtomImpl.AromaticSubset.Boron);
                    g.markDelocalised();
                    break;
                case 'c':
                    addAtom(AtomImpl.AromaticSubset.Carbon);
                    g.markDelocalised();
                    break;
                case 'n':
                    addAtom(AtomImpl.AromaticSubset.Nitrogen);
                    g.markDelocalised();
                    break;
                case 'o':
                    addAtom(AtomImpl.AromaticSubset.Oxygen);
                    g.markDelocalised();
                    break;
                case 'p':
                    addAtom(AtomImpl.AromaticSubset.Phosphorus);
                    g.markDelocalised();
                    break;
                case 's':
                    addAtom(AtomImpl.AromaticSubset.Sulfur);
                    g.markDelocalised();
                    break;


                // D/T for hydrogen isotopes - non-standard but OpenSMILES spec
                // says it's possible. The D and T here are automatic converted
                // to [2H] and [3H].
                case 'H':
                    if (strict)
                        throw new InvalidSmilesException("hydrogens should be specified in square brackets - '[H]'",
                                                         buffer);
                    addAtom(AtomImpl.EXPLICIT_HYDROGEN);
                    break;
                case 'D':
                    if (strict)
                        throw new InvalidSmilesException("deuterium should be specified as a hydrogen isotope - '[2H]'",
                                                         buffer);
                    addAtom(AtomImpl.DEUTERIUM);
                    break;
                case 'T':
                    if (strict)
                        throw new InvalidSmilesException("tritium should be specified as a hydrogen isotope - '[3H]'",
                                                         buffer);
                    addAtom(AtomImpl.TRITIUM);
                    break;

                // bracket atom
                case '[':
                    addAtom(readBracketAtom(buffer));
                    break;

                // ring bonds
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    ring(c - '0', buffer);
                    break;
                case '%':
                    int num = buffer.getNumber(2);
                    if (num < 0)
                        throw new InvalidSmilesException("a number (<digit>+) must follow '%':", buffer);
                    ring(num, buffer);
                    break;

                // bond/dot
                case '-':
                    bond = Bond.SINGLE;
                    break;
                case '=':
                    bond = Bond.DOUBLE;
                    break;
                case '#':
                    bond = Bond.TRIPLE;
                    break;
                case '$':
                    bond = Bond.QUADRUPLE;
                    break;
                case ':':
                    g.markDelocalised();
                    bond = Bond.AROMATIC;
                    break;
                case '/':
                    bond = Bond.UP;
                    break;
                case '\\':
                    bond = Bond.DOWN;
                    break;
                case '.':
                    bond = Bond.DOT;
                    break;

                // branching
                case '(':
                    if (stack.empty())
                        throw new InvalidSmilesException("cannot open branch - there were no previous atoms:",
                                                         buffer);
                    stack.push(stack.peek());
                    break;
                case ')':
                    if (stack.size() < 2)
                        throw new InvalidSmilesException("closing of an unopened branch:",
                                                         buffer);
                    stack.pop();
                    break;

                // termination
                case '\t':
                case ' ':
                case '\n':
                case '\r':
                    return;

                default:
                    throw new InvalidSmilesException("unexpected character:", buffer);
            }
        }
    }


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
