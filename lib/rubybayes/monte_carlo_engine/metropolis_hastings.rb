module Rubybayes
  module MonteCarloEngine
  
    module MetropolisHastings_Log
      def rho(y)
        f(y) - f(@xt) + g(@xt) - g(y) 
      end
      def eval(y)
        Math.log(rand) < rho(y)
      end	 
    end
    
    module MetropolisHastings_Normal
      def rho(y)   
        (f(y) / f(@xt)) * (g(@xt) / g(y))        
      end      
      def eval(y)
        Kernel.rand < rho(y)
      end 
    end
    
    module MetropolisHastings_RandomWalk
      def g(x)   
        1.0
      end      
      def sample_from_g
       @g.rng(@xt)
      end 
    end
    
    module MetropolisHastings_Independent
      def g(x)   
        @g.pdf(x)
      end      
      def sample_from_g
        @g.rng()
      end 
    end
    
    class MetropolisHastings            
     
      def initialize(args)
        @f = args[:f] # function(s): pdf
        @g = args[:g] # function(s): rng, pdf
        @xt = args[:start]
        
        if args.fetch(:random_walk, false) == true
          extend Rubybayes::MonteCarloEngine::MetropolisHastings_RandomWalk
        else
          extend Rubybayes::MonteCarloEngine::MetropolisHastings_Independent
        end
                
        if args.fetch(:log, false) == true 
          extend Rubybayes::MonteCarloEngine::MetropolisHastings_Log
        else
          extend Rubybayes::MonteCarloEngine::MetropolisHastings_Normal
        end
        
        @accept = 0
        @counter = 0
      end

      def accepted
        @accept.to_f / @counter.to_f 
      end
      
      def f(x)
        @f.pdf(x)
      end
      
      def sample
        @counter += 1
        y = sample_from_g
        if (eval(y)) then
          @accept += 1
          @xt = y
        end
        @xt
      end	
      
    end

  end
end