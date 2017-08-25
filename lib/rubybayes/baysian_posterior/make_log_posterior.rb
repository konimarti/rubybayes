require "rubybayes/baysian_posterior/minimization"

module Rubybayes

  module MakeLogPosterior
  
    def proposal_density
      raise "You should implement this"
    end
  
    def log_likelihood(x)
      raise "You should implement this"
    end
  
    def log_priors(x)
      raise "You should implement this"
    end
  
    def starting_point
      raise "You should implement this"
    end
    
    def pdf(x)
      log_likelihood(x) + log_priors(x)
    end
  
    def one_parameter_mode(lower=-1000.0, upper=1000.0)
        #Brent for one-parameter cases
        min = Rubybayes::Minimization::Brent.new(lower,upper, proc {|x| -self.pdf(x)})
        min.iterate
        min.x_minimum	    
    end
  
    def mode
      #Nelder-Mead for multi-parameter cases
      min = Rubybayes::Minimization::NelderMead.new(proc {|x| -self.pdf(x)}, starting_point)
      i = 0
      while min.converging? 
        i += 1
        min.iterate 
      end   
      min.x_minimum       
    end
  
    def generate_engine
      Rubybayes::MonteCarloEngine::MetropolisHastings.new(
        f: self,
        g: proposal_density,        
        log: true,
        random_walk: random_walk,
        start: starting_point
      )
    end
    
    def engine
      @engine ||= generate_engine
    end
    
    def sample
      engine.sample
    end
  
    def random_walk
      false
    end
  
  end
  
  #'prepend' this module
  module EnableRandomWalk
    def random_walk
      true
    end
  end
  
end