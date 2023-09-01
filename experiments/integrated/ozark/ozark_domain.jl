# Domain setup
# For soil column
nelements = 40
zmin = FT(-2)
zmax = FT(0)
println("depth of soil domain is  "*string(zmin))
println("dx of soil domain is  "*string(round((zmax-zmin)./nelements,digits=3)))
land_domain =
    LSMSingleColumnDomain(; zlim = (zmin, zmax), nelements = nelements)

# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
n_stem = Int64(1);
n_leaf = Int64(1);
println("numbrer of stem and leaf compartments are  "*string(n_stem)*" "*string(n_leaf))
h_stem = FT(9) # m, from Wang et al.
println("height of stem compartment is "*string(h_stem))
h_leaf = FT(9.5) # m from Wang et al.
println("height of leaf compartment is "*string(h_leaf))
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2];
compartment_surfaces = [zmax, h_stem, h_stem + h_leaf];
# For soil column
