module Rubybayes
  module MonteCarloEngine
  
    module MetropolisHastingsLog
    
      def rho(y)
        fy = f(y)
        fx = f(@xt)
        gy = g(y)
        gx = g(@xt)
        (fy - fx + gx - gy) 
      end
      
      def eval(y)
        Math.log(rand) < rho(y)
      end	 
      
    end
    
    module MetropolisHastingsNormal
    
      def rho(y)
        fy = f(y)
        fx = f(@xt)
        gy = g(y)
        gx = g(@xt)       
        (fy / fx) * (gx / gy)        
      end
      
      def eval(y)
        Kernel.rand < rho(y)
      end	   
      
    end
    
    class MetropolisHastings            
      attr_reader :accept 
      def initialize(args)
        @f = args[:f] # function(s): pdf
        @g = args[:g] # function(s): rng, pdf
        @xt = args[:start]
        
        @random_walk = args.fetch(:random_walk, false)
                
        if args.fetch(:log, false) 
          extend Rubybayes::MonteCarloEngine::MetropolisHastingsLog
        else
          extend Rubybayes::MonteCarloEngine::MetropolisHastingsNormal
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
      
      def g(x)
        if @random_walk
          1.0
        else
          @g.pdf(x)
        end
      end

      def sample_from_g
        if @random_walk
          @g.rng(@xt)
        else
          @g.rng
        end
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