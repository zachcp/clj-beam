; Try to outline a Clojure-based version of the Beam Smiles Parsing Library
; https://github.com/johnmay/beam
;

(ns clj-beam.core
  (:gen-class)
  (:require
   [schema.core :as s]
   [plumbing.core :as p :include-macros true]
   [plumbing.fnk.pfnk :as pfnk :include-macros true]
   [plumbing.graph :as graph :include-macros true]
   [plumbing.map :as map])
  ;(:import [uk.ac.ebi.beam Atom Graph Parser CharBuffer])
  )


;; Graph g = Graph.fromSmiles("CCO");
;; String smi = g.toSmiles();
;; Graph Calls Parser
;; Parser calls CharBuffer


;; basic Goal is to get a semidecnt parser to get a datastructure

(defn roundtrip [smiles]
  (let [s   (Graph/fromSmiles smiles)
        ret (.toSmiles s)]
    ret))

;; try to model my calculations after Prismatic using Schem to describe datashapes and plumbing
;; to describe each set of transformations that the parsing undergoes.

(defn parse-smiles [smiles]
    (seq smiles))

(parse-smiles "CCO")

;https://github.com/johnmay/beam/blob/d6ab8c8c916e51ecbc715a9e7b537f85431dc05d/core/src/main/java/uk/ac/ebi/beam/Parser.java
(defn readSmiles [s] )

; Use the same Grpah Model as Beam for holding information about the Chemical Structure
(defrecord cAtom  isotope element aromatic charge hydrogens atomClass subset)
(defrecord cGraph atoms edges topologies order size)
(defrecord cParser stack graph ringbonds arrangement configurations)

(def charmap {  '*' :Unknown
                ;note that disambiguation of the list is needed
                'B' :Bromine
                'B' :Boron
                'C' :Chlorine
                'C' :Carbon
                'N' :Nitrogen
                'O' :Oxygen
                'P' :Phosphorus
                'S' :Sulfur
                'F' :Fluorine
                'I' :Iodine
                ;note that aromatinc is not taken care of
                'b' :Boron
                'c' :Carbon
                'n' :Nitrogen
                'o' :Oxygen
                'p' :Phosphorus
                's' :Sulfur
                'D' :DEUTERIUM
                'T' :TRITIUM
                ;note buffer needs dismiguation
                '\[' :BUFFER
                '0' :RING
                '1' :RING
                '2' :RING
                '3' :RING
                '4' :RING
                '5' :RING
                '6' :RING
                '7' :RING
                '8' :RING
                '9' :RING
                '\%' :temp
                   ; BUFFER.getnumber(2)
                    ;case '%':
                    ;int num = buffer.getNumber(2);
                    ;if (num < 0)
                    ;    throw new InvalidSmilesException("number (<digit>+) must follow '%'", buffer);
                    ;ring(num);
                    ;break;

                ;// bond/dot
                '\-' :Bond.SINGLE;
                '\=' :Bond.DOUBLE;
                '#' :Bond.TRIPLE;
                '$' :Bond.QUADRUPLE;
                ':' :Bond.AROMATIC;
                '\/' :Bond.UP
                '\\\\' :Bond.DOWN;
                '.' :Bond.DOT;
                ;// branching
                '\(' :xx
                    ;if (stack.empty()) throw new InvalidSmilesException("cannot open branch - no previous atom",
                    ;                                     buffer);
                    ;stack.push(stack.peek());
                    ;break;
                '\)' :cvcv
                    ;if (stack.size() < 2)
                    ;    throw new InvalidSmilesException("Closing of an unopened branch",

                ;// termination
                '\t' :end
                ' '  :end
                '\n' :end
                '\r' :end })



(for [meth (.getMethods Parser)
      :let [name (.getName meth)]]
  name)

(println .getMethods CharBuffer)
