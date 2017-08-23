require "test_helper"
require "rubystats"

class Posterior
	def initialize
		@n = 9 
		@w = 6
		@prior = Rubystats::NormalDistribution.new(0.7,1.0)
	end
	
	def pdf(x)
		Math.log(likelihood(x)) + Math.log(prior(x))
	end
	
	def likelihood(x)
		Rubystats::BinomialDistribution.new(@n, x).pdf(@w)
	end
	
	def prior(x)
		@prior.pdf(x)
	end
end

class BetaDistWrapper
	def initialize(a,b)
		@beta = Rubystats::BetaDistribution.new(a,b)
	end
	def method_missing(name,*args,&block)
		@beta.send(name, *args, &block)
	end
  def pdf(x)
    Math.log(@beta.pdf(x))
  end
	def rng(x=nil)
		@beta.icdf(rand)
	end
end

class SimpleBaysianInferenceTest < Minitest::Test
  def test_simple_test
    #create monte carlo engine for sampling a markov chain 
    engine = Rubybayes::MonteCarloEngine::MetropolisHastings.new(
      f: Posterior.new,
      g: BetaDistWrapper.new(6,3),
      start: 0.5,
      log: true,
      random_walk: true
    )
    #run monte carlo simulation
    experiment = Rubybayes::MonteCarloSimulation::Simulation.new do
      burn_in 1000
      iterations 5000
      sample { engine.sample }  
    end
    #perform analysis on markov chain
    m = Rubybayes::MonteCarloSimulation.extract_measurements(experiment.run.get_chains)            
    
    assert_in_epsilon 0.6363, m[0].mean, 0.1   
  end
end
