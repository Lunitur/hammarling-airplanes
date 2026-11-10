{
  description = "Julia i LaTeX razvojna okruženja za Hammarlingov algoritam";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { nixpkgs, ... }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "aarch64-darwin" ];
    in
    {
      devShells = nixpkgs.lib.genAttrs systems (system:
        let
          pkgs = import nixpkgs { inherit system; };
          jupyter = pkgs.python3.withPackages (ps: [ ps.notebook ]);
          tex = pkgs.texliveBasic.withPackages (ps: with ps; [
            amsmath
            amsfonts
            babel-croatian
            booktabs
            cm-super
            enumitem
            geometry
            hyphen-croatian
            latexmk
          ]);
        in
        {
          default = pkgs.mkShell {
            packages = [ pkgs.julia-bin jupyter ];
            JULIA_PROJECT = "@.";
            JUPYTER = "${jupyter}/bin/jupyter";
            IJULIA_NODEFAULTKERNEL = "1";
            shellHook = ''
              julia --project=@. -e '
                using Pkg
                Pkg.instantiate()
                project = dirname(Base.active_project())
                Pkg.activate("hammarling-ijulia-v$(VERSION.major).$(VERSION.minor)"; shared=true)
                if !haskey(Pkg.project().dependencies, "IJulia")
                    Pkg.add("IJulia")
                else
                    Pkg.instantiate()
                end
                using IJulia
                IJulia.installkernel("Julia Hammarling", "--project=$project";
                    specname="julia-hammarling",
                    env=Dict("JULIA_LOAD_PATH" => "@:$(dirname(Base.active_project())):@stdlib"))
              '
            '';
          };

          latex = pkgs.mkShell {
            packages = [ tex pkgs.texlab ];
          };
        });
    };
}
