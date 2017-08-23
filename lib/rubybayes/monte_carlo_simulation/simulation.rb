module Rubybayes

  module MonteCarloSimulation
  
    class Simulation
    
      attr_accessor :result
      
      def initialize(&block)
        @results = {}
        #set defaults
        @burn_in  = 0
        @iterations = 1000
        @chains = 1
        @calculate = Proc.new {|x| x}
        @engine = nil
        #set user input
        self.instance_eval(&block)
      end
      
      # public functions
      
      def run        
        @chains.times do |nchain| 
          @results[nchain]=[]
          run_chain(nchain)
        end
        self
      end      
      
      def get_merged_chains
        merge_chains(@results)
      end      
            
      def get_chains
        @results        
      end
      
      private
      
      # DSL
      
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
            
      # run simulation  
         
      def run_chain(nchain)        
        raise ArgumentError.new("'sample' not set for Monte Carlo simulation") if @sample.nil? 
        @burn_in.times {|i| @sample.call() } if @burn_in > 0        
        @iterations.times {|i| @results[nchain] << @calculate.call( @sample.call() )}
      end     
          
    end 
    
    # module functions    
    
    def self.merge_chains(results)
      ret = []
      results.each_value {|vals| ret.concat(vals) }
      {0 => ret}    
    end         
     
    def self.run(&block)
      Rubybayes::MonteCarloSimulation::Simulation.new(&block).run.get_chains
    end
    
    def self.run_and_combine(&block)
      Rubybayes::MonteCarloSimulation::Simulation.new(&block).run.get_merged_chains
    end    
    
  end
  
end