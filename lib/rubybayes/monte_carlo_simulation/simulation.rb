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
        @burn_in = x.to_i
      end
      
      def iterations(x)
        @iterations = x.to_i
      end
      
      def chains(x)
        @chains = x.to_i
      end
      
      def generate_sampler(&block)
        @engine = block
      end
      
      def sample(&block)
        @sample = block
      end
      
      def calculate(&block)
        @calculate = block
      end
         
      def run_chain
        raise "not implemented yet"
      end
         
      def run
        raise ArgumentError.new("'sample' not set for Monte Carlo simulation") if @sample.nil? 
        @burn_in.times {|i| @sample.call() } if @burn_in > 0        
        @iterations.times {|i| @results << @calculate.call( @sample.call() )}
        @results
      end   
          
    end
    
    def self.run(&block)
      Rubybayes::MonteCarloSimulation::Simulation.new(&block).run
    end
    
  end
  
end