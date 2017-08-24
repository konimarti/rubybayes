require "test_helper"
require "rubystats"

class Posterior
	def initialize
		@n = 9 
		@w = 6
		@prior = Rubystats::BetaDistribution.new(6,3)
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


class FullyFeaturedPosterior
  include Rubybayes::MakeLogPosterior
  
  def initialize
  	@n = 9 
		@w = 6
		@prior = Rubystats::BetaDistribution.new(6,3)
  end
  
  def proposal_density
    BetaDistWrapper.new(6,3)
  end
  
  def log_likelihood(x)
    Math.log(Rubystats::BinomialDistribution.new(@n, x).pdf(@w))
  end
  
  def log_priors(x)
    Math.log(@prior.pdf(x))
  end
  
  def starting_point
    0.5
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
      random_walk: false
    )
    
    #run monte carlo simulation
    experiment = Rubybayes::MonteCarloSimulation::Simulation.new do
      burn_in 1000
      iterations 10000
      sample { engine.sample }  
    end
    
    #perform analysis on markov chain
    m = Rubybayes::MonteCarloSimulation.extract_measurements(experiment.run)            
    
    assert_in_epsilon 0.6667, m[0].mean, 0.01   
  end
  
  
  def test_make_log_posterior
    posterior = FullyFeaturedPosterior.new
    mode = posterior.one_parameter_mode(0.0,1.0)
    
    1000.times{posterior.sample}
    ret = []
    10000.times{ret << posterior.sample}
    
    m = Rubybayes::MonteCarloSimulation::Measurement.new(ret)  
    
    #puts "mode = #{mode}"
    assert_in_epsilon 0.6875, mode, 0.001  
    
    #puts "average = #{m.mean}"    
    assert_in_epsilon 0.6667, m.mean, 0.01   
  end
  
  
end
