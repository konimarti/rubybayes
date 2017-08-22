module Rubybayes

  module MonteCarloSimulation
  
    class Simulation
    
      attr_accessor :result
      
      def initialize(&block)
        @results = []
        #set defaults
        @burn_in  = 0
        @iterations = 10000
        @calculate = Proc.new {|x| x}
        #set user input
        self.instance_eval(&block)
      end
      
      def burn_in(x)
        @burn_in = x
      end
      
      def iterations(x)
        @iterations = x
      end
      
      def sample(&block)
        @sample = block
      end
      
      def calculate(&block)
        @calculate = block
      end
            
      def run
        raise ArgumentError.new("'sample' not set for Monte Carlo simulation") if @sample.nil? 
        @burn_in.times {|i| @sample.call() }        
        @iterations.times {|i| @results << @calculate.call( @sample.call() )}
        @results
      end   
          
    end
    
    def self.run(&block)
      Rubybayes::MonteCarloSimulation::Simulation.new(&block).run
    end
    
  end
  
end