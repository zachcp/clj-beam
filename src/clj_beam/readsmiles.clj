;; Smile Specification
;; http://www.daylight.com/dayhtml/doc/theory/theory.smiles.htmlsmiles.clj

; Try to outline a Clojure-based version of the Beam Smiles Parsing Library
; https://github.com/johnmay/beam
;
;
;
;
; Basic Idea:
;   1. Get a Sequence of characters from a string
;   2. Use that Sequence to define elements
;   3. Use that Sequence to define bonds
;   4. Use that Sequence to Define Rings
;   5. Use That Sequence to define sterochemistry
;


(ns clj-beam.readsmiles
  (:gen-class)
  (:require
   [clojure.string :as string]
   [schema.core :as s]
   [plumbing.core :as p :include-macros true]
   [plumbing.fnk.pfnk :as pfnk :include-macros true]
   [plumbing.graph :as graph :include-macros true]
   [plumbing.map :as map]
   [clj-beam.schemas :as schemas]))



;;; Master Function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(def readsmiles
  "A graph specifying the parsing of a smiles string"
  {
   :sqs  (p/fnk [smi] (string->keys smi))
   :atms (p/fnk [sqs] (into [] (map symbol->Atom sqs)))

   ; need:
   ; 1: find isotopes
   ; 2: find charges
   ; 2: find find branches
   ; 3: allow bonding between branch and moelcuel before it
   ; 4: some global state to avoid doubling up

   :sbnds (p/fnk [sqs] (detect-singlebonds sqs))
   :dbnds (p/fnk [sqs] (detect-doublebonds sqs))
   :tbnds (p/fnk [sqs] (detect-triplebonds sqs))
   :bnds (p/fnk [sqs] (str "needs to be done"))
   :g    (p/fnk [atms sbonds dbonds tbonds]
                {:atoms  [atms]
                 :bonds  [ stms sbonds dbonds tbonds ]})})



;; Funcitons
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn string->keys [smi]
  "Take a Smiles String and get a sequence of parsed smiles units"
  (let [conv (fn [x] (cond
                       ( = (apply str x) "Br") :Br
                       ( = (apply str x) "Cl") :Cl
                       :else (keyword (str (first x)))))]
        (->> smi
             (partition 2 1 '(:Padding)); partition into characters and use :Padding
             (filter #(not= (first %) \r)) ;BR abd Cl are special cases
             (filter #(not= (first %) \l))
             (map conv)
             (vec))))


(defn symbol->Atom
  "Map Symbol to Atom Creation"
  [charvect]
  (for [ idx (range (count v)) ]
    (let [ sym (get charvect idx)]
       (cond
        (in? [:Br :B :Cl :C :N :O :P :S :F :I] sym )
             {:element ( sym elementmap) :aromatic :No :index idx}
        (in? [:c :n :o :p :s ] sym )
             {:element (sym elementmap) :aromatic :Yes}
        (in? [:H :D :T] sym )
             {:element (sym elementmap) :aromatic :No}
        :else "There is a problem" ))))


(defn detect-singlebonds
  "Parse Sequence and Generate Sequence Bonds. Note: This will NOT currently work on anyhting inside of brackets"
   [charvect]
   (let [ ;keep track of adjacent indeces
          part  (partition 2 1 '(:Padding) charvect)
          idx   (partition 2 1 '(:Padding) (range (count charvect)))
          inter (interleave part idx)]
     (for [ [p i] (partition 2 inter ) ]
        (cond
           (and (in? (keys elementmap) (first p))
                (in? (keys elementmap) (second p)))
             {:atom1 (first i) :atom2 (second i) :order :Single :aromatic :No} ;return a map corresponding to the bond
                ))))

(defn detect-rings
  "Detect Ring Indices by Location of the Numbers: Note Need to check for aromatic or not!!!"
   [charvect]
     (let [i     (range 15)  ;assuming arbitrary maximum number of rings
           locs  (fn [x] (get-indices (keyword (str x)) charvect))
           locs1 (map locs i)
           ;note there couls be More than 2 if there are many nubmers
           bonds (filter #(= (count %) 2) locs1)]

         (for [b bonds]
             {:atom1 (- (first b) 1) :atom2 (- (second b) 1)}  )))

(defn detect-doublebonds
  "Detect Double Bonds: Note that this method doe snot see things in brackets"
   [charvect]
     (let [dbs (get-indices := charvect)]
         (for [d dbs]
             {:atom1 (- d 1) :atom2 (+ d 1) :order :Double}  )))


(defn detect-triplebonds
  "Detect Triple Bonds: Note that this method doe snot see things in brackets"
   [charvect]
     (let [dbs (get-indices :# charvect)]
         (for [d dbs]
             {:atom1 (- d 1) :atom2 (+ d 1) :order :Double}  )))


;; Scratch
(string->keys smi3)
sequences2
(detect-single-bonds sequences2)
(detect-rings sequences2)
(detect-doublebonds sequences2)


;; QuickChecks
;; Check the Basic Parsing
(def smiles (string/split  (slurp "data/smiles.txt") #"\n"))
(def smi1 (first smiles))
(def smi2 (second smiles))
(def smi3 (nth smiles 3))
(def smi6 (nth smiles 6))
smi3
(= (string->keys smi1) '(:C :C))
(= (string->keys smi2) '(:O := :C := :O))
(= (string->keys smi3) '(:C :C :N :( :C :C :) :C :C))
(= (string->keys (nth smiles 11)) '(:F :/ :C := :C :\ :F))

(def sequences (string->keys smi3))
(def sequences2 (string->keys smi6))
(def atoms (symbol->Atom sequences))



;;; Helper Functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def elementmap {
  :* :Unknown    :Br :Bromine   :B  :Boron
  :Cl :Chlorine  :C :Carbon     :N :Nitrogen
  :O :Oxygen     :P :Phosphorus :S :Sulfur
  :F :Fluorine   :I :Iodine     :b :Boron
  :c :Carbon     :n :Nitrogen   :o :Oxygen
  :p :Phosphorus :s :Sulfur     :H :Hydrogen
  :D :DEUTERIUM  :T :TRITIUM })

(defn in?
  "true if seq contains elm"
  [seq elm]
  (some #(= elm %) seq))

(defn notin?
  "tru if seq does not contian element"
  [seq elm]
  (not (in? seq elm)))

(defn get-indices [x coll]
  (map first
    (filter #(= (second %) x)
    (map-indexed vector coll))))

