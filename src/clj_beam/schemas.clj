(ns clj-beam.schemas
  (:require
   [schema.core :as s]
   [plumbing.core :as p :include-macros true]
   [plumbing.fnk.pfnk :as pfnk :include-macros true]
   [plumbing.graph :as graph :include-macros true]
   [plumbing.map :as map])
  )


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
