(defproject dynamic-notebook/functional-numerics "0.1.1"
  :description "Functional Numerical Methods in Clojure(script)."
  :url "http://github.com/dynamic-notebook/functional-numerics"
  :scm {:name "git" :url "http://github.com/dynamic-notebook/functional-numerics"}
  :dependencies [[org.clojure/clojure "1.10.1" :scope "provided"]
                 [org.clojure/clojurescript "1.10.773" :scope "provided"]
                 [org.clojure/core.match "1.0.0"]]
  :profiles {:dev {:plugins [[lein-cljsbuild "1.1.8"
                              :exclusions [org.clojure/clojurescript]]
                             [lein-doo "0.1.11"]]
                   :cloverage {:ns-exclude-regex [#"sicmutils.rules"
                                                  #"sicmutils.simplify"]}
                   :repl-options {:nrepl-middleware
                                  [cider.piggieback/wrap-cljs-repl]}
                   :dependencies [[org.clojure/test.check "1.1.0"]
                                  [same/ish "0.1.4"]
                                  [cider/piggieback "0.5.0"]
                                  [lein-doo "0.1.11"]]}
             :test {:jvm-opts ["-Xmx512m"]
                    :dependencies [[org.clojure/test.check "1.0.0"]
                                   [criterium "0.4.5"]]}}
  :aliases {"test-cljs" ["doo" "node" "test" "once"]}
  :cljsbuild {:builds
              {:test
               {:source-paths ["src" "test"]
                :compiler
                {:main sicmutils.runner
                 :optimizations :none
                 :target :nodejs
                 :output-dir "target/main"
                 :output-to "target/main/main.js"}}}})
