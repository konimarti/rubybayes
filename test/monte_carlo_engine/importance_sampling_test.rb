require "test_helper"
require "rubystats"

class ShiftedExponential
  def initialize
    @exp = Rubystats::ExponentialDistribution.new(2)
  end
  def pdf(x)
    @exp.pdf(x-3)
  end
  def rng
    @exp.rng + 3
  end
end

class ImportanceSamplingTest < Minitest::Test
  def test_gaussian_tail_probability
  
    #create monte carlo engine for sampling  
    engine = Rubybayes::MonteCarloEngine::ImportanceSampling.new(
      f: Rubystats::NormalDistribution.new(0.0,1.0),
      #g: ShiftedExponential.new
      g: Rubystats::NormalDistribution.new(4.0,1.0)
    )
    
    #run monte carlo simulation for P(Y > 3) with Y ~ N(0,1)
    experiment = Rubybayes::MonteCarloSimulation::Simulation.new do
      iterations 100000
      sample { engine.sample }
      calculate {|x| h=x[0]; w=x[1]; ((h>3.0)?(1.0):(0.0))*w}      
    end
    #perform analysis on markov chain
    m = Rubybayes::MonteCarloSimulation.extract_measurements(experiment.run)            
    
    puts "integral = #{m[0].mean}"
    
    assert_in_delta 0.001349, m[0].mean, 0.0001   
  end
end
