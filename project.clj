(defproject clj-beam "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [prismatic/plumbing "0.3.1"]
                 [prismatic/schema "0.2.2"]
                 [uk.ac.ebi.beam/beam-core "0.6"]
                 [uk.ac.ebi.beam/beam-func "0.6"]
                 ;[uk.ac.ebi.beam/beam-core "0.7-SNAPSHOT"]
                 ;[uk.ac.ebi.beam/beam-func "0.7-SNAPSHOT"]
                 ]
  :main ^:skip-aot clj-beam.core
  :target-path "target/%s"
  :repositories [["ebi-repo" "http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo" ]
                 ["ebi-repo-snapshots" "http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo" ]]
  :profiles {:uberjar {:aot :all}})
