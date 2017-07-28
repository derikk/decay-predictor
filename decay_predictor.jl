# Predicts nuclear binding energies and decays of nuclides
# Copyright © 2017 by Derik Kauffman. Licensed under MIT license.

# To identify elements by their names
element_symbols = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]
element_names = ["Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminium", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon", "Caesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium", "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinum", "Gold", "Mercury", "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", "Francium", "Radium", "Actinium", "Thorium", "Protactinium", "Uranium", "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium", "Einsteinium", "Fermium", "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium", "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium", "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium", "Flerovium", "Moscovium", "Livermorium", "Tennessine", "Oganesson"]

# Masses in MeV/c²
const mass = Dict(:p => 938.272081, :n => 939.565413, :e => 0.510999, :α => 3727.379378)

# Compiled sets of coefficients
averaged_coefficients = (15.78, 18., 0.712, 23.5, 11.8, 0.5)
tinker_coefficients = (15.78, 18.03, 0.713, 23.5, 11.8, 0.5)
wiki_ls1_coefficients = (15.8, 18.3, 0.714, 23.2, 12, 0.5)
wiki_ls2_coefficients = (15.76, 17.81, 0.711, 23.702, 34, 0.75)
wiki_rohlf_coefficients = (15.75, 17.8, 0.711, 23.7, 11.18, 0.5)
tipler_coefficients = (15.76, 17.23, 0.75, 23.2, 12., 0.5)  # aA was written as 93.2 but that's way off
vahid_coefficients = (15.519, 17.476, 0.674, 24.576, 12, 0.5)
hyperphys_coefficients = (15.75, 17.8, 0.711, 23.7, 11.18, 0.5)

aV, aS, aC, aA, aP, kP = averaged_coefficients





# I present to you… the semi-empirical mass formula! Courtesy of Carl Friedrich von Weizsäcker.
# Check it out: https://en.wikipedia.org/wiki/Semi-empirical_mass_formula
binding_E(A, Z) = aV*A - aS*A^(2/3) - aC*Z^2*A^(-1/3) - aA*(A-2Z)^2/A + δ(A, Z) + magic(A, Z)
binding_per_nucleon(A, Z) = binding_E(A, Z)/A

# Bonus points if Z and/or N are even.
δ(A, Z) = -((A+1) % 2) * (Z%2 * 2 - 1) * aP/A^kP

# Yes, this function is magic. See https://en.wikipedia.org/wiki/Magic_number_(physics)
function magic(A, Z)
  magic_nums = [2, 8, 20, 28, 50, 82]  # Add more?
  aM = 0.8  # Need to find a good value for this.
  aM*((Z in magic_nums) + (A-Z in magic_nums))  # Divide by √A?
end

Qval_α_decay(A, Z) = (2mass[:p] + 2mass[:n] - binding_E(A, Z)) - (mass[:α] - binding_E(A-4, Z-2))
Qval_β⁻_decay(A, Z) = (mass[:n] - binding_E(A, Z)) - (mass[:p] + mass[:e] - binding_E(A, Z+1))
Qval_β⁺_decay(A, Z) = (mass[:p] - binding_E(A, Z)) - (mass[:n] + mass[:e] - binding_E(A, Z-1))

println(binding_E(50, 26))  # Actual:  417.70
println(binding_E(235, 92)) # Actual: 1783.87
println(binding_E(238, 92)) # Actual: 1801.69

println(Qval_α_decay(235, 92)) # Actual: 4.68
println(Qval_α_decay(238, 92)) # Actual: 4.27
println(Qval_β⁻_decay(62, 26)) # Actual: 2.53
println(Qval_β⁺_decay(53, 26)) # Actual: 2.72

function decay_modes(A, Z)
  α_decays = Qval_α_decay(A, Z) > 0
  β⁻_decays = Qval_β⁻_decay(A, Z) > 0
  β⁺_decays = Qval_β⁺_decay(A, Z) > 0
  e_captures = Qval_β⁺_decay(A, Z) > -2mass[:e] # Instead of emitting positron, capture electron

  print(element_names[Z], "-", A)
  if α_decays || β⁻_decays || e_captures
    if α_decays
      print(" will α-decay to form ", element_names[Z-2], "-", A-4)
      if β⁻_decays
        print(" or β⁻-decay to form ", element_names[Z+1], "-", A)
      elseif e_captures
        if β⁺_decays
          print(" or undergo electron capture or β⁺-decay to form ", element_names[Z-1], "-", A)
        else
          print(" or decay by electron capture to form ", element_names[Z-1], "-", A)
        end
      end
    elseif β⁻_decays
      print(" will β⁻-decay to form ", element_names[Z+1], "-", A)
    elseif β⁺_decays
      print(" will undergo electron capture or β⁺-decay to form ", element_names[Z-1], "-", A)
    else
      print(" will decay by electron capture to form ", element_names[Z-1], "-", A)
    end
  else
    print(" is stable")
  end
  println(".")
end


#=== Testing ===#

decay_modes(238, 92)

for A in 52:60
  decay_modes(A, 26)
end
println(Qval_β⁻_decay(60, 26))  # Not stable! Half life of a couple million years.
#=for A in 79:86
  decay_modes(A, 36)
end=#

for A in 202:209
  decay_modes(A, 82)
  println("Pb-", A, " BE: ", binding_E(A, 82))
end
println("Pb-208 BE: ", binding_E(208, 82)) # 1636.43 MeV
println("Pb-208 Qval: ", Qval_β⁻_decay(208, 82)) # Stable
for A in 39:48
  decay_modes(A, 20)
end
println(Qval_β⁺_decay(41, 20)) # EC, 0.421 MeV
println(binding_E(41, 20)) # 350.41 MeV
println(Qval_β⁻_decay(46, 20)) # Stable
println("Ca-48 BE: ", binding_E(48, 20)) #
println("Ca-48 Qval: ", Qval_β⁻_decay(48, 20)) # Practically stable
#=for A in 56:65
  decay_modes(A, 28)
end=#
decay_modes(132, 55)
println("Cs-132 BE: ", binding_E(132, 55)) #
println("Cs-132 Qval-: ", Qval_β⁻_decay(132, 55))
println("Cs-132 Qval+: ", Qval_β⁺_decay(132, 55))
