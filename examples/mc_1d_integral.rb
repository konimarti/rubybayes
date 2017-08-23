require "rubybayes"

include Rubybayes::MonteCarloSimulation

#define function to evaluate
def f(x)
  -1.0 * x**2 + 1
end

#define simulation
experiment = Simulation.new do
  iterations 50000
  chains 4
  sample { 2.0 * (Kernel.rand - 0.5) }
  calculate {|x| f(x)}
end

#run Monte Carlo Simulation
results = experiment.run

#merge chains and create Measurement object
measurement = Measurement.new( Rubybayes::MonteCarloSimulation.merge_chains(results)[0] )

#output
integral = 2.0 * measurement.mean

puts "Calculated integral from -1.0 to 1.0 = #{integral}"
puts "Correct integral = 1.33333"