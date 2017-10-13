with import <nixpkgs> {};

let
   hinge = callPackage ./default.nix {};
in

stdenv.mkDerivation {
  name = "hinge-env";
  buildInputs = [hinge];
}

