{
  gcc,
  stdenv,
  cmake,
  hdf5,
  python3Packages,
  callPackage,
  pkgs,
}:

let
   python = callPackage ./requirements.nix { inherit pkgs; };
in

  stdenv.mkDerivation {
    name = "hinge";

    src = ./.;

    buildInputs = with python3Packages; [

      gcc
      stdenv
      cmake
      hdf5

      (boost.override { enableStatic = true;}) 

      #numpy
      scipy

      python.packages.networkx
      python.packages.easydev
      python.packages.colormap
      python.packages.ujson
    ];

    #   checkPhase = ''
    #   export NO_TESTS_OVER_WIRE=1
    #   export PYTHONDONTWRITEBYTECODE=1

    #   #flake8 src/
    #   #py.test --cov=src -cov-report term-missing
    #   #coverage html
    # '';

    # setting up the right compiler for boost
    # BOOST_ROOT = "${boost}";
    # BOOST_LIBRARYDIR = "${boost}";
    # CMAKE_C_COMPILER = "${gcc49}";
    # CMAKE_CXX_COMPILER = "${gcc49}";

    # shellHook =
    # ''
    # rm -rf build
    # source $stdenv/setup
    # export BOOST_INCLUDEDIR=${boost.dev}/include
    # export BOOST_LIBRARYDIR=${boost.out}/lib
    # mkdir build && cd build
    # cmake .. -DCMAKE_INSTALL_PREFIX=../inst -DBOOST_LIBRARYDIR=${boost.out}/lib
    # # bash ./utils/build.sh
    # '';
}
