# generated using pypi2nix tool (version: 1.8.0)
# See more at: https://github.com/garbas/pypi2nix
#
# COMMAND:
#   pypi2nix -V 3 -e ujson -e colormap -e easydev -e networkx
#

{ pkgs ? import <nixpkgs> {}
}:

let

  inherit (pkgs) makeWrapper;
  inherit (pkgs.stdenv.lib) fix' extends inNixShell;

  pythonPackages =
  import "${toString pkgs.path}/pkgs/top-level/python-packages.nix" {
    inherit pkgs;
    inherit (pkgs) stdenv;
    python = pkgs.python3;
  };

  commonBuildInputs = [];
  commonDoCheck = false;

  withPackages = pkgs':
    let
      pkgs = builtins.removeAttrs pkgs' ["__unfix__"];
      interpreter = pythonPackages.buildPythonPackage {
        name = "python3-interpreter";
        buildInputs = [ makeWrapper ] ++ (builtins.attrValues pkgs);
        buildCommand = ''
          mkdir -p $out/bin
          ln -s ${pythonPackages.python.interpreter}               $out/bin/${pythonPackages.python.executable}
          for dep in ${builtins.concatStringsSep " "               (builtins.attrValues pkgs)}; do
            if [ -d "$dep/bin" ]; then
              for prog in "$dep/bin/"*; do
                if [ -f $prog ]; then
                  ln -s $prog $out/bin/`basename $prog`
                fi
              done
            fi
          done
          for prog in "$out/bin/"*; do
            wrapProgram "$prog" --prefix PYTHONPATH : "$PYTHONPATH"
          done
          pushd $out/bin
          ln -s ${pythonPackages.python.executable} python
          popd
        '';
        passthru.interpreter = pythonPackages.python;
      };
    in {
      __old = pythonPackages;
      inherit interpreter;
      mkDerivation = pythonPackages.buildPythonPackage;
      packages = pkgs;
      overrideDerivation = drv: f:
        pythonPackages.buildPythonPackage (drv.drvAttrs // f drv.drvAttrs);
      withPackages = pkgs'':
        withPackages (pkgs // pkgs'');
    };

  python = withPackages {};

  generated = self: {

    "colorama" = python.mkDerivation {
      name = "colorama-0.3.9";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/e6/76/257b53926889e2835355d74fec73d82662100135293e17d382e2b74d1669/colorama-0.3.9.tar.gz"; sha256 = "48eb22f4f8461b1df5734a074b57042430fb06e1d61bd1e11b078c0fe6d7a1f1"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [ ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = licenses.bsdOriginal;
        description = "Cross-platform colored terminal text.";
      };
    };



    "colormap" = python.mkDerivation {
      name = "colormap-1.0.1";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/68/8f/9cf6bd10f0d2c4522504412ace0faeaac67f8e6280af6965a88f825c712d/colormap-1.0.1.tar.gz"; sha256 = "d361296c796d0868e0fbed2506b0a774f739819d2446a2ce7f8dd4e4c228477a"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [ ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = licenses.lgpl2;
        description = "Utilities to ease manipulation of matplotlib colormaps and color codecs (e.g., hex2rgb)";
      };
    };



    "decorator" = python.mkDerivation {
      name = "decorator-4.1.2";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/bb/e0/f6e41e9091e130bf16d4437dabbac3993908e4d6485ecbc985ef1352db94/decorator-4.1.2.tar.gz"; sha256 = "7cb64d38cb8002971710c8899fbdfb859a23a364b7c99dab19d1f719c2ba16b5"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [ ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = licenses.bsdOriginal;
        description = "Better living through Python with decorators";
      };
    };



    "easydev" = python.mkDerivation {
      name = "easydev-0.9.35";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/03/41/4838f6c6bd89ba5fe67973e9d1f20a64786c98e6779e10a9ec1adb1ec7b5/easydev-0.9.35.tar.gz"; sha256 = "8663a7012c34b6be9f9ceb2c39cee7725b6516a95db41083ae34cad8a756950a"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [
      self."colorama"
      self."pexpect"
    ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = "License :: OSI Approved";
        description = "Common utilities to ease the development of Python packages";
      };
    };



    "networkx" = python.mkDerivation {
      name = "networkx-2.0";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/73/58/0add7d81cf64958f7a062aa287237364eb0a0959bf7a817f708d5c25d043/networkx-2.0.zip"; sha256 = "cd5ff8f75d92c79237f067e2f0876824645d37f017cfffa5b7c9678cae1454aa"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [
      self."decorator"
    ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = licenses.bsdOriginal;
        description = "Python package for creating and manipulating graphs and networks";
      };
    };



    "pexpect" = python.mkDerivation {
      name = "pexpect-4.2.1";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/e8/13/d0b0599099d6cd23663043a2a0bb7c61e58c6ba359b2656e6fb000ef5b98/pexpect-4.2.1.tar.gz"; sha256 = "3d132465a75b57aa818341c6521392a06cc660feb3988d7f1074f39bd23c9a92"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [
      self."ptyprocess"
    ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = licenses.isc;
        description = "Pexpect allows easy control of interactive console applications.";
      };
    };



    "ptyprocess" = python.mkDerivation {
      name = "ptyprocess-0.5.2";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/51/83/5d07dc35534640b06f9d9f1a1d2bc2513fb9cc7595a1b0e28ae5477056ce/ptyprocess-0.5.2.tar.gz"; sha256 = "e64193f0047ad603b71f202332ab5527c5e52aa7c8b609704fc28c0dc20c4365"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [ ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = "";
        description = "Run a subprocess in a pseudo terminal";
      };
    };



    "ujson" = python.mkDerivation {
      name = "ujson-1.35";
      src = pkgs.fetchurl { url = "https://pypi.python.org/packages/16/c4/79f3409bc710559015464e5f49b9879430d8f87498ecdc335899732e5377/ujson-1.35.tar.gz"; sha256 = "f66073e5506e91d204ab0c614a148d5aa938bdbf104751be66f8ad7a222f5f86"; };
      doCheck = commonDoCheck;
      buildInputs = commonBuildInputs;
      propagatedBuildInputs = [ ];
      meta = with pkgs.stdenv.lib; {
        homepage = "";
        license = licenses.bsdOriginal;
        description = "Ultra fast JSON encoder and decoder for Python";
      };
    };

  };
  overrides = import ./requirements_override.nix { inherit pkgs python; };
  commonOverrides = [

  ];

in python.withPackages
   (fix' (pkgs.lib.fold
            extends
            generated
            ([overrides] ++ commonOverrides)
         )
   )