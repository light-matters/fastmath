(ns build
  "Compiles fastmath's Java sources into target/classes.

  Used as the `:deps/prep-lib` step when fastmath is a tools.deps git
  dependency (consumers run `clj -X:deps prep`), and runnable directly with
  `clj -T:build compile-java`.

  LBFGSBJava is a git submodule, and tools.deps git dependencies never fetch
  submodules. When its sources are absent, the pinned submodule commit is
  cloned (https) into target/LBFGSBJava and compiled from there."
  (:require
   [clojure.java.io :as io]
   [clojure.tools.build.api :as b]))

(def class-dir "target/classes")

(def lbfgsb
  {:url "https://github.com/genmeblog/LBFGSBJava.git"
   :sha "be2fdae9550f459930b2d075c7b477fedefadf80"   ; = `git submodule status` LBFGSBJava
   :checkout "target/LBFGSBJava"})

(defn- sh! [& args]
  (let [{:keys [exit]} (b/process {:command-args (vec args)})]
    (when-not (zero? exit)
      (throw (ex-info (str "command failed: " (pr-str args)) {:exit exit})))))

(defn- lbfgsb-src
  "Source directory of LBFGSBJava: the submodule if checked out, otherwise
  a clone of the pinned commit."
  []
  (if (.exists (io/file "LBFGSBJava" "src"))
    "LBFGSBJava/src"
    (let [{:keys [url sha checkout]} lbfgsb]
      (when-not (.exists (io/file checkout "src"))
        (b/delete {:path checkout})
        (sh! "git" "clone" "--quiet" url checkout)
        (sh! "git" "-C" checkout "checkout" "--quiet" sha))
      (str checkout "/src"))))

(defn compile-java [_]
  (b/javac {:src-dirs ["src" (lbfgsb-src)]
            :class-dir class-dir
            :basis (b/create-basis {:project "deps.edn"})
            :javac-opts ["--release" "8" "-Xlint:unchecked"]}))
