module Rubybayes
  module MonteCarloEngine
  
    class ImportanceSampling            
     
      def initialize(args)
        @f = args[:f] # function(s): pdf
        @g = args[:g] # function(s): rng, pdf        
      end

      def f(x)
        @f.pdf(x)
      end
      
      def g(x)
        @g.pdf(x)
      end
          
      def sample_from_g
        @g.rng
      end
          
      def sample
        xt = sample_from_g
		    [xt, f(xt)/g(xt)]
      end	
      
    end

  end
end