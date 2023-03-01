require_relative "lib/partitioned_array/lib/line_db"


require 'bigdecimal'
require 'matrix'
#require 'matrix/dsl'


class Particle
  attr_accessor :position, :mass, :splitting_potential_power, :wave_function

  def initialize(position, mass, splitting_potential_power = 1, wave_function = Quantum::Qubit.new)
    @position = position
    @mass = mass
    @splitting_potential_power = splitting_potential_power
    if wave_function.is_a?(Matrix)
      @wave_function = Quantum::Qubit.new(wave_function[0, 0], wave_function[1, 0])
    elsif wave_function.is_a?(Quantum::Qubit)
      @wave_function = wave_function
    else
      raise ArgumentError, "Invalid wave function type: #{wave_function.class}"
    end
  end
  def kinetic_energy(velocity, mass = @mass)
    BigDecimal("0.5") * mass * velocity ** 2
  end

  def to_str
    Matrix.column_vector([@position, @mass, @splitting_potential_power])
  end

  def tensor_product(other_particle)
    new_position = @position.to_f * other_particle.position.to_f
    new_mass = @mass.to_f * other_particle.mass.to_f
    new_splitting_potential_power = @splitting_potential_power.to_f + other_particle.splitting_potential_power.to_f
    new_wave_function = @wave_function.tensor_product(other_particle.wave_function)
    Particle.new(new_position, new_mass, new_splitting_potential_power, new_wave_function[0,0], new_wave_function[1,0])
  end
end

class NeutriParticle
  attr_accessor :mass, :splitting_potential_power, :position, :wave_function, :qubit


  def initialize(mass = BigDecimal(0), splitting_potential_power = BigDecimal(2), position = BigDecimal(0))
    @mass = mass
    @splitting_potential_power = splitting_potential_power
    @position = position
    @qubit = Quantum::Qubit.new
    @qubit.alpha = Complex(1, 0)
    @qubit.beta = Complex(0, 0)
    @state = @qubit.to_vector
    @wave_function = @state # initialize wave_function with state
  end
  
  
  def tensor_product(other_particle)
    new_qubit = @qubit.tensor_product(other_particle.qubit)
    new_position = @position.kronecker_product(other_particle.position)
    NeutriParticle.new(@mass * other_particle.mass, @splitting_potential_power + other_particle.splitting_potential_power, new_position).tap do |new_particle|
      new_particle.qubit = new_qubit
      new_particle.wave_function = new_qubit.to_vector
    end
  end
  def probability(state)
    # Calculate the probability of measuring the particle in the given state
    (state.conjugate * @qubit.to_vector).abs2
  end

  def superposition(other_particle)
    # Create a new qubit that represents the superposition of the two particles
    new_qubit = Quantum::Qubit.new
    new_qubit.alpha = @qubit.alpha * other_particle.qubit.alpha
    new_qubit.beta = @qubit.beta * other_particle.qubit.beta
    
    new_qubit_other = Quantum::Qubit.new
    new_qubit_other.alpha = @qubit.alpha * other_particle.qubit.beta
    new_qubit_other.beta = @qubit.beta * other_particle.qubit.alpha
    
    # Update the qubits of the particles
    @qubit = new_qubit
    other_particle.qubit = new_qubit_other
    
    # Update the wave function of the particles based on the new qubit state
    @wave_function = new_qubit.to_vector
    other_particle.wave_function = new_qubit_other.to_vector
  end


  def wave_length
    if @mass == 0
      0
    else
      BigDecimal('6.62607004e-34') / BigDecimal('1.60217662e-19') / BigDecimal(299792458) / BigDecimal(@mass.abs)
    end
  end

  

  def momentum
    if @mass == 0
      0
    else
      BigDecimal('6.62607004e-34') / BigDecimal('1.60217662e-19') / wave_length
    end
  end

  def spin
    @mass == 0 ? BigDecimal(0) : @mass / @splitting_potential_power
  end

  def kinetic_energy(velocity, mass = @mass)
    BigDecimal("0.5") * mass * velocity ** 2
  end

  def impact(z)
    if (@mass.abs < BigDecimal("1e-10")) && (@splitting_potential_power.abs < BigDecimal("1e-10")) || ((@mass == 0) && (@splitting_potential_power == 0)) || z <= 0
      1e-10 # changed from "(0)"
    else
      (@mass / (z + @splitting_potential_power)) * z
    end
  end

  def additive_factors(m, c)
    c == 0 ? 1 : (BigDecimal(m) / BigDecimal(c * c)).ceil
  end

  class ResonanceFrequency
    def initialize(neutri_particle)
      @neutri_particle = neutri_particle
    end

    def calculate
      @neutri_particle.mass / @neutri_particle.splitting_potential_power
    end
  end

  class Entropy
    def initialize(neutri_particle)
      @neutri_particle = neutri_particle
    end

    def calculate
      additive_factor = @neutri_particle.additive_factors(@neutri_particle.mass, @neutri_particle.splitting_potential_power)
      BigDecimal(@neutri_particle.mass * @neutri_particle.splitting_potential_power * @neutri_particle.impact(BigDecimal(1)) * additive_factor)
    end
  end

  class Entanglement
    def initialize(neutri_particle)
      @neutri_particle = neutri_particle
    end

    def calculate(other_particle)
      additive_factor = @neutri_particle.additive_factors(@neutri_particle.spin, @neutri_particle.mass)
      @neutri_particle.spin * other_particle.position * additive_factor
    end
  end
end

class NeutriParticleInteraction
  attr_accessor :particles, :distance

  def initialize(*particles, distance: BigDecimal(1))
    @particles = particles
    @distance = distance
  end

  def probability_of_interaction
    energy_factor = BigDecimal(0)
    mass_factor = BigDecimal(1)
    spin_factor = BigDecimal(1)
  
    @particles.each do |particle|
      if particle.is_a?(NeutriParticle)
        energy_factor += particle.kinetic_energy(BigDecimal(1)) 
        spin_factor *= particle.spin
      else
        mass_factor *= particle.mass
      end
    end
  
    additive_factor = BigDecimal(1)
    if mass_factor > BigDecimal(0) && spin_factor > BigDecimal(0)
      spin_additive_factor = (spin_factor / (mass_factor ** BigDecimal("0.3333333333333333333333333333333333333333333333333333333333333333333"))).ceil
      mass_additive_factor = (mass_factor / BigDecimal("1e3")).ceil
      ke_additive_factor = (energy_factor / BigDecimal("1e3")).ceil
      additive_factor = spin_additive_factor * mass_additive_factor * ke_additive_factor
    end
  
    dist_factor = BigDecimal(1) / (BigDecimal(@distance ** 2))
    probability = dist_factor * energy_factor / (mass_factor + energy_factor) * additive_factor
    return probability
  end
end

class ProbabilityOfInteraction
  def initialize(interaction)
    @interaction = interaction
  end

  def calculate
    return @interaction.probability_of_interaction
  end
end

module Quantum
  class Qubit
    attr_accessor :alpha, :beta

    def initialize(alpha = Complex(1, 0), beta = Complex(0, 0))
      @alpha = alpha
      @beta = beta
    end

    def measure
      rand_num = rand()
      if rand_num <= @alpha.abs2
        @alpha = Complex(1, 0)
        @beta = Complex(0, 0)
        return 0
      else
        @alpha = Complex(0, 0)
        @beta = Complex(1, 0)
        return 1
      end
    end
    
    def tensor_product(other_qubit)
      # compute the tensor product between this qubit and another qubit
      alpha_beta = [@alpha*@beta, other_qubit.alpha*other_qubit.beta]
      beta_alpha = [@beta*other_qubit.alpha, @alpha*other_qubit.beta]
      Quantum::Qubit.new(alpha_beta.sum, beta_alpha.sum)
    end

    def to_matrix
      Matrix.column_vector([[@alpha],[@beta]])
    end

    def to_vector
      Matrix.column_vector([@alpha, @beta])
    end

    def to_s
      "alpha: #{@alpha}, beta: #{@beta}"
    end

    def hadamard_gate
      # apply the Hadamard gate to the qubit
      @alpha, @beta = (1 / Math.sqrt(2)) * (@alpha + @beta), (1 / Math.sqrt(2)) * (@alpha - @beta)
    end

    def pauli_x_gate
      # apply the Pauli X gate to the qubit
      @alpha, @beta = @beta, @alpha
    end

    def pauli_y_gate
      # apply the Pauli Y gate to the qubit
      @alpha, @beta = Complex(-@beta.imag, @alpha.imag), Complex(@beta.real, -@alpha.real)
    end

    def pauli_z_gate
      # apply the Pauli Z gate to the qubit
      @beta = -@beta
    end
  end
end






# Create two particles
particle1 = Particle.new(BigDecimal(1), BigDecimal(1))
particle2 = Particle.new(BigDecimal(2), BigDecimal(2))

# Create two neutrino particles
neutri_particle1 = NeutriParticle.new(BigDecimal(3), BigDecimal(3))
neutri_particle2 = NeutriParticle.new(BigDecimal(4), BigDecimal(4))

# Calculate resonance frequency, entropy, and entanglement for each neutrino particle
rf1 = NeutriParticle::ResonanceFrequency.new(neutri_particle1).calculate
rf2 = NeutriParticle::ResonanceFrequency.new(neutri_particle2).calculate

entropy1 = NeutriParticle::Entropy.new(neutri_particle1).calculate
entropy2 = NeutriParticle::Entropy.new(neutri_particle2).calculate

entanglement1 = NeutriParticle::Entanglement.new(neutri_particle1).calculate(particle1)
entanglement2 = NeutriParticle::Entanglement.new(neutri_particle2).calculate(particle2)

# Create two neutrino interactions with distance of 3 units
interaction1 = NeutriParticleInteraction.new(neutri_particle1, particle1, neutri_particle2, particle2, distance: BigDecimal(3))
interaction2 = NeutriParticleInteraction.new(neutri_particle2, particle2, neutri_particle1, particle1, distance: BigDecimal(3))

# Calculate probability of interaction for each interaction
prob1 = interaction1.probability_of_interaction
prob2 = interaction2.probability_of_interaction

# Print out all the results
puts "Particle 1 position: #{particle1.position}, mass: #{particle1.mass}"
puts "Particle 2 position: #{particle2.position}, mass: #{particle2.mass}"

puts "Neutrino particle 1 mass: #{neutri_particle1.mass}, splitting potential power: #{neutri_particle1.splitting_potential_power}"
puts "Neutrino particle 2 mass: #{neutri_particle2.mass}, splitting potential power: #{neutri_particle2.splitting_potential_power}"

puts "Neutrino particle 1 resonance frequency: #{rf1}"
puts "Neutrino particle 2 resonance frequency: #{rf2}"

puts "Neutrino particle 1 entropy: #{entropy1}"
puts "Neutrino particle 2 entropy: #{entropy2}"

puts "Neutrino particle 1 entanglement with particle 1: #{entanglement1}"
puts "Neutrino particle 2 entanglement with particle 2: #{entanglement2}"

puts "Neutrino interaction 1 probability of interaction: #{prob1}"
puts "Neutrino interaction 2 probability of interaction: #{prob2}"

# Instantiate Neutri particles
neutri1 = NeutriParticle.new(BigDecimal("1"), BigDecimal("2"))
neutri2 = NeutriParticle.new(BigDecimal("1"), BigDecimal("2"))
neutri3 = NeutriParticle.new(BigDecimal("1"), BigDecimal("2"))
neutri4 = NeutriParticle.new(BigDecimal("1"), BigDecimal("2"))

# Calculate probability of interaction for all pairs of particles with the same distance
distance = BigDecimal("2")
interactions = []
[neutri1, neutri2, neutri3, neutri4].combination(2).each do |particle1, particle2|
  interaction = NeutriParticleInteraction.new(particle1, particle2, distance: distance)
  probability = ProbabilityOfInteraction.new(interaction).calculate
  interactions << [particle1, particle2, probability]
end

# Print out the results
puts "Interactions with distance #{distance}:"
interactions.each do |particle1, particle2, probability|
  puts "#{particle1.inspect} and #{particle2.inspect} have a probability of interaction of #{probability}"
end

neutrino1 = NeutriParticle.new(BigDecimal(1), BigDecimal(2))
neutrino2 = NeutriParticle.new(BigDecimal(2), BigDecimal(3))

interaction = NeutriParticleInteraction.new(neutrino1, neutrino2, distance: BigDecimal(1))
probability = ProbabilityOfInteraction.new(interaction).calculate

puts "Probability of interaction: #{probability}"

# create NeutriParticle and Particle instances
neutri_particle_1 = NeutriParticle.new(BigDecimal('1e-20'), BigDecimal('1e-20'))
neutri_particle_2 = NeutriParticle.new(BigDecimal('1e-20'), BigDecimal('1e-20'))
particle_1 = Particle.new(BigDecimal('1e-18'), BigDecimal('1e-25'))
particle_2 = Particle.new(BigDecimal('1e-18'), BigDecimal('1e-25'))

# create NeutriParticleInteraction and ProbabilityOfInteraction instances
interaction_1 = NeutriParticleInteraction.new(neutri_particle_1, neutri_particle_2, particle_1, particle_2, distance: 0)
interaction_2 = NeutriParticleInteraction.new(neutri_particle_1, neutri_particle_2, particle_1, particle_2, distance: BigDecimal('1e-25'))
interaction_3 = NeutriParticleInteraction.new(neutri_particle_1, neutri_particle_2, particle_1, particle_2, distance: BigDecimal('1e-30'))
probability_1 = ProbabilityOfInteraction.new(interaction_1)
probability_2 = ProbabilityOfInteraction.new(interaction_2)
probability_3 = ProbabilityOfInteraction.new(interaction_3)

# calculate and print probability of interaction for each instance
puts "Probability of interaction for interaction_1: #{probability_1.calculate}"
puts "Probability of interaction for interaction_2: #{probability_2.calculate}"
puts "Probability of interaction for interaction_3: #{probability_3.calculate}"

neutron = Particle.new(BigDecimal("1.00866491600"), BigDecimal("1.675e-27"))
neutri = NeutriParticle.new(BigDecimal("1e-24"), BigDecimal("1e-10"))

interaction = NeutriParticleInteraction.new(neutron, neutri, distance: 2)
probability = ProbabilityOfInteraction.new(interaction).calculate

puts "The probability of interaction between a neutron and a neutrino at a distance of 2 is #{probability}."

# Calculating Resonance Frequency
puts "Resonance Frequency"
neutri_particle_a = NeutriParticle.new(BigDecimal("1e-6"), BigDecimal(2))
resonance_frequency = NeutriParticle::ResonanceFrequency.new(neutri_particle_a).calculate
puts "Resonance frequency of NeutriParticle with mass #{neutri_particle_a.mass} and splitting potential power #{neutri_particle_a.splitting_potential_power}: #{resonance_frequency}"

# Calculating Entropy
puts "Entropy"
neutri_particle_b = NeutriParticle.new(BigDecimal("1e-6"), BigDecimal(2))
entropy = NeutriParticle::Entropy.new(neutri_particle_b).calculate
puts "Entropy of NeutriParticle with mass #{neutri_particle_b.mass} and splitting potential power #{neutri_particle_b.splitting_potential_power}: #{entropy}"

# Calculating Entanglement
puts "Entanglement"
neutri_particle1 = NeutriParticle.new(BigDecimal("1e-6"), BigDecimal(2))
neutri_particle2 = NeutriParticle.new(BigDecimal("1e-5"), BigDecimal(3))
entanglement = NeutriParticle::Entanglement.new(neutri_particle1).calculate(neutri_particle2)
puts "Entanglement between NeutriParticle1 with mass #{neutri_particle1.mass} and splitting potential power #{neutri_particle1.splitting_potential_power} and NeutriParticle2 with mass #{neutri_particle2.mass} and splitting potential power #{neutri_particle2.splitting_potential_power}: #{entanglement}"

# Calculating Probability of Interaction
puts "Probability of Interaction"
neutri_particle1 = NeutriParticle.new(BigDecimal("1e-6"), BigDecimal(2))
neutri_particle2 = NeutriParticle.new(BigDecimal("1e-5"), BigDecimal(3))
particle1 = Particle.new(BigDecimal("1e-6"), BigDecimal(2))
particle2 = Particle.new(BigDecimal("1e-5"), BigDecimal(3))
interaction = NeutriParticleInteraction.new(neutri_particle1, neutri_particle2, particle1, particle2, distance: 2)
probability = ProbabilityOfInteraction.new(interaction).calculate
puts "Probability of interaction between NeutriParticle1 with mass #{neutri_particle1.mass} and splitting potential power #{neutri_particle1.splitting_potential_power} and NeutriParticle2 with mass #{neutri_particle2.mass} and Particle1 with mass #{particle1.mass} and Particle2 with mass #{particle2.mass} and splitting potential power #{particle2.splitting_potential_power} at a distance of #{interaction.distance}: #{probability}"

# Define the particles
neutri1 = NeutriParticle.new(BigDecimal("1e-9"), BigDecimal("1e-9"))
neutri2 = NeutriParticle.new(BigDecimal("1e-9"), BigDecimal("1e-9"))
particle1 = Particle.new(BigDecimal("1"), BigDecimal("1"))
particle2 = Particle.new(BigDecimal("1"), BigDecimal("1"))

# Define the range of distances to check
min_distance = BigDecimal("1e-10")
max_distance = BigDecimal("1e-5")
step_size = BigDecimal("1e-9")

# Iterate through different distances and calculate the probability of interaction
max_probability = 0
max_distance = 5
current_distance = min_distance
while current_distance <= max_distance
  interaction = NeutriParticleInteraction.new(neutri1, neutri2, particle1, particle2, distance: current_distance)
  probability = ProbabilityOfInteraction.new(interaction).calculate
  if probability > max_probability
    max_probability = probability
    max_distance = current_distance
  end
  current_distance += step_size
end

# Output the distance and probability with the highest value
puts "Maximum probability of interaction: #{max_probability}"
puts "Distance with maximum probability: #{max_distance}"

particle1 = Particle.new(BigDecimal("0.5"), BigDecimal("100"))
particle2 = Particle.new(BigDecimal("1"), BigDecimal("100"))

neutri1 = NeutriParticle.new(BigDecimal("0.01"), BigDecimal("0.5"))
neutri2 = NeutriParticle.new(BigDecimal("0.01"), BigDecimal("0.5"))

interaction1 = NeutriParticleInteraction.new(neutri1, neutri2, particle1, particle2, distance: BigDecimal("0.1"))
interaction2 = NeutriParticleInteraction.new(neutri1, particle1, distance: BigDecimal("0.01"))

probabilities = []
[interaction1, interaction2].each do |interaction|
  prob = ProbabilityOfInteraction.new(interaction).calculate
  probabilities << { interaction: interaction, probability: prob }
end

# sort by probability in descending order
probabilities.sort_by! { |p| -p[:probability] }

# print information for highest probability interaction
highest_probability = probabilities.first
puts "Highest probability interaction:"
puts "Particles: #{highest_probability[:interaction].particles.map(&:class)}"
puts "Distance: #{highest_probability[:interaction].distance}"
puts "Probability: #{highest_probability[:probability]}"

# print information for second highest probability interaction
second_highest_probability = probabilities[1]
puts "Second highest probability interaction:"
puts "Particles: #{second_highest_probability[:interaction].particles.map(&:class)}"
puts "Distance: #{second_highest_probability[:interaction].distance}"
puts "Probability: #{second_highest_probability[:probability]}"


# Define particles
proton = Particle.new(BigDecimal('0.00000167'), BigDecimal('1.00728'), BigDecimal('0'))
neutri1 = NeutriParticle.new(BigDecimal('0.0000000001'), BigDecimal('2'), BigDecimal('0'))
neutri2 = NeutriParticle.new(BigDecimal('0.0000000001'), BigDecimal('2'), BigDecimal('0'))

# Define interaction
interaction = NeutriParticleInteraction.new(proton, neutri1, neutri2, distance: BigDecimal('1'))

# Find highest probability
max_probability = BigDecimal('0')
max_distance = nil
max_particles = []
(1..10).each do |dist|
  interaction.distance = BigDecimal(distance)
  probability = ProbabilityOfInteraction.new(interaction).calculate
  if probability > max_probability
    max_probability = probability
    max_distance = dist
    max_particles = [proton] + interaction.particles.select{|p| p.is_a?(NeutriParticle)}
  end
end

# Output highest probability and distance
puts "Maximum probability: #{max_probability}"
puts "At distance: #{max_distance}"
puts "Particle details:"
max_particles.each do |particle|
  puts "Mass: #{particle.mass}, Splitting Potential Power: #{particle.splitting_potential_power}, Position: #{particle.position}"
end

n = NeutriParticle.new(BigDecimal("3"),BigDecimal(2))
impact_notation = n.impact(BigDecimal("999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999"))
puts "Impact parameter notation: #{impact_notation}"



# Example 1
n_particle = NeutriParticle.new(BigDecimal("0.1"), BigDecimal("2"), BigDecimal("1"))
n_interaction = NeutriParticleInteraction.new(n_particle, n_particle, distance: BigDecimal("3"))
p = ProbabilityOfInteraction.new(n_interaction)
puts p.calculate # Output: 4.4e-8

# Example 2
n_particle1 = NeutriParticle.new(BigDecimal("0.3"), BigDecimal("1"), BigDecimal("0.1"))
n_particle2 = NeutriParticle.new(BigDecimal("0.5"), BigDecimal("3"), BigDecimal("2"))
n_interaction = NeutriParticleInteraction.new(n_particle1, n_particle2, distance: BigDecimal("5"))
p = ProbabilityOfInteraction.new(n_interaction)
puts p.calculate # Output: 2.61e-10

# Example 3
n_particle1 = NeutriParticle.new(BigDecimal("0.3"), BigDecimal("1"), BigDecimal("0.1"))
n_particle2 = Particle.new(BigDecimal("1"), BigDecimal("0.1"), BigDecimal("0"))
n_interaction = NeutriParticleInteraction.new(n_particle1, n_particle2, distance: BigDecimal("2"))
p = ProbabilityOfInteraction.new(n_interaction)
puts p.calculate # Output: 5.5e-3

# Example 4
n_particle1 = NeutriParticle.new(BigDecimal("0"), BigDecimal("2"), BigDecimal("0"))
n_particle2 = NeutriParticle.new(BigDecimal("1"), BigDecimal("2"), BigDecimal("1"))
n_interaction = NeutriParticleInteraction.new(n_particle1, n_particle2, distance: BigDecimal("1"))
p = ProbabilityOfInteraction.new(n_interaction)
puts p.calculate # Output: (0)

# Example 5
n_particle1 = NeutriParticle.new(BigDecimal("0"), BigDecimal("0"), BigDecimal("0"))
n_particle2 = NeutriParticle.new(BigDecimal("1"), BigDecimal("1"), BigDecimal("0"))
n_interaction = NeutriParticleInteraction.new(n_particle1, n_particle2, distance: BigDecimal("2"))
p = ProbabilityOfInteraction.new(n_interaction)
puts p.calculate # Output: (0)

# Example 6
n_particle1 = NeutriParticle.new(BigDecimal("0"), BigDecimal("0"), BigDecimal("1"))
n_particle2 = Particle.new(BigDecimal("1"), BigDecimal("1"), BigDecimal("1"))
n_interaction = NeutriParticleInteraction.new(n_particle1, n_particle2, distance: BigDecimal("5"))
p = ProbabilityOfInteraction.new(n_interaction)
puts p.calculate # Output: (0)



neutri_particle1 = NeutriParticle.new(BigDecimal("1"), BigDecimal("2"), BigDecimal("0.5"))
neutri_particle2 = NeutriParticle.new(BigDecimal("0.5"), BigDecimal("1"), BigDecimal("1.5"))
particle = Particle.new(BigDecimal("3"), BigDecimal("2"), BigDecimal("1"))
interaction = NeutriParticleInteraction.new(neutri_particle1, neutri_particle2, particle, distance: BigDecimal("1"))
probability = ProbabilityOfInteraction.new(interaction).calculate
p probability #=> 0.0793650793...

interaction = NeutriParticleInteraction.new(neutri_particle1, neutri_particle2, particle, distance: BigDecimal("2"))
probability = ProbabilityOfInteraction.new(interaction).calculate
p probability #=> 0.0198412698...

interaction = NeutriParticleInteraction.new(neutri_particle1, neutri_particle2, particle, distance: BigDecimal("3"))
probability = ProbabilityOfInteraction.new(interaction).calculate
p probability #=> 0.0082644628...


##########

# Define a particle with non-zero mass and splitting potential power
particle1 = NeutriParticle.new(BigDecimal("1"), BigDecimal("2"))

# Define a second particle with zero mass and splitting potential power
particle2 = NeutriParticle.new(BigDecimal("0"), BigDecimal("0"))

# Define a third particle with non-zero mass and splitting potential power
particle3 = NeutriParticle.new(BigDecimal("2"), BigDecimal("1"))

# Define a distance of 2 between each of the particles
distance = BigDecimal("2")

# Create an interaction between all three particles at the given distance
interaction = NeutriParticleInteraction.new(particle1, particle2, particle3, distance: distance)

# Calculate the probability of interaction between all three particles
probability = ProbabilityOfInteraction.new(interaction).calculate

puts "The probability of interaction between particles 1, 2, and 3 at a distance of #{distance} is #{probability}."


# Example of wave-particle duality
particle1 = NeutriParticle.new(BigDecimal(1))
particle2 = NeutriParticle.new(BigDecimal(2))
interaction = NeutriParticleInteraction.new(particle1, particle2, distance: BigDecimal('2'))
probability_calculator = ProbabilityOfInteraction.new(interaction)

puts "Wave length of particle 1: #{particle1.wave_length}"
puts "Wave length of particle 2: #{particle2.wave_length}"
puts "Momentum of particle 1: #{particle1.momentum}"
puts "Momentum of particle 2: #{particle2.momentum}"
puts "Probability of interaction: #{probability_calculator.calculate}"


require 'bigdecimal'

# define a NeutriParticle with a very large mass and very small splitting potential power
np = NeutriParticle.new(BigDecimal("1e20"), BigDecimal("1e-10"))

# calculate the wavelength and momentum of the particle
wave_length = np.wave_length
momentum = np.momentum

puts "Wave length: #{wave_length}"
puts "Momentum: #{momentum}"



####




##########

# Example of wave-particle duality
particle1 = NeutriParticle.new(BigDecimal(1), 2)
particle2 = NeutriParticle.new(BigDecimal(2), 2)
particle3 = Particle.new(BigDecimal(1), 2)
interaction = NeutriParticleInteraction.new(particle1, particle2, particle3, distance: BigDecimal('2'))
probability_calculator = ProbabilityOfInteraction.new(interaction)

puts "Wave length of particle 1: #{particle1.wave_length}"
puts "Wave length of particle 2: #{particle2.wave_length}"
puts "Momentum of particle 1: #{particle1.momentum}"
puts "Momentum of particle 2: #{particle2.momentum}"
puts "Probability of interaction: #{probability_calculator.calculate}"
puts "Distance: #{interaction.distance}"
puts "Particle 1 position: #{particle1.position}"
puts "Particle 2 position: #{particle2.position}"
puts "Particle 3 position: #{particle3.position}"
puts "Particle 1 mass: #{particle1.mass}"
puts "Particle 2 mass: #{particle2.mass}"
puts "Particle 3 mass: #{particle3.mass}"
puts "Particle 1 splitting potential power: #{particle1.splitting_potential_power}"
puts "Particle 2 splitting potential power: #{particle2.splitting_potential_power}"
puts "Particle 3 splitting potential power: #{particle3.splitting_potential_power}"
puts particle1.impact(42)
puts particle2.impact(42)
#puts particle3.impact(42)


puts "####################"

neutrino = NeutriParticle.new(0.05, 0.1)
target = Particle.new(BigDecimal(1000), BigDecimal(100))
impact = neutrino.impact(target.position)
probability = ProbabilityOfInteraction.new(NeutriParticleInteraction.new(neutrino, target)).calculate
puts "The impact of the neutrino on the target is #{impact}."
puts "The probability of the neutrino interacting with the target is #{probability}."


# Define the particles
particle1 = NeutriParticle.new(BigDecimal("0.1"), BigDecimal("1"))
particle2 = NeutriParticle.new(BigDecimal("0.2"), BigDecimal("2"))
particle3 = NeutriParticle.new(BigDecimal("0.3"), BigDecimal("3"))

# Create an interaction between the particles
interaction = NeutriParticleInteraction.new(particle1, particle2, particle3, distance: BigDecimal("2"))

# Calculate the probability of interaction
prob = ProbabilityOfInteraction.new(interaction).calculate
puts "Probability of interaction: #{prob}"

# Print the particles' state before superposition
puts "Particle 1 state before superposition: #{particle1.wave_function}"
puts "Particle 2 state before superposition: #{particle2.wave_function}"
puts "Particle 3 state before superposition: #{particle3.wave_function}"

# Create a superposition between the particles
particle1.superposition(particle2)
particle1.superposition(particle3)

# Print the particles' state after superposition
puts "Particle 1 state after superposition: #{particle1.wave_function}"
puts "Particle 2 state after superposition: #{particle2.wave_function}"
puts "Particle 3 state after superposition: #{particle3.wave_function}"

puts "####################"

particle1 = NeutriParticle.new(BigDecimal(1), BigDecimal(1))
particle2 = NeutriParticle.new(BigDecimal(1), BigDecimal(1))
particle3 = NeutriParticle.new(BigDecimal(1), BigDecimal(1))

puts "Particle 1 state before superposition: #{particle1.wave_function}"
puts "Particle 2 state before superposition: #{particle2.wave_function}"
puts "Particle 3 state before superposition: #{particle3.wave_function}"

particle1.superposition(particle2)
particle2.superposition(particle3)

puts "\nAfter superposition:\n"

puts "Particle 1 state after superposition: #{particle1.wave_function}"
puts "Particle 2 state after superposition: #{particle2.wave_function}"
puts "Particle 3 state after superposition: #{particle3.wave_function}"

puts "####################"

qubits = []

# create 64 qubits and put them into superposition
1024.times do
  qubit = Quantum::Qubit.new
  qubit.hadamard_gate
  qubits << qubit
end

# print the superposition states of the qubits
qubits.each_with_index do |qubit, index|
  puts "Qubit #{index + 1} state: #{qubit.to_vector}"
end
puts "####################"

# Create 64 qubits in the |0> state
qubits = Array.new(5) { Quantum::Qubit.new }

# Apply the Hadamard gate to each qubit
qubits.each(&:hadamard_gate)

# Compute the tensor product of all qubits to get the 64-qubit superposition state
superposition = qubits.reduce { |result, qubit| result.tensor_product(qubit) }

puts superposition.alpha
puts superposition.beta

puts "####################"
particles = []

20000.times do
  mass = BigDecimal(rand(1..999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999))
  splitting_potential_power = BigDecimal(rand(1..999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999))
  position = BigDecimal(rand(1..99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999))
  particles << NeutriParticle.new(mass, splitting_potential_power, position)
end

interaction = NeutriParticleInteraction.new(*particles)

probability = ProbabilityOfInteraction.new(interaction).calculate
p probability
