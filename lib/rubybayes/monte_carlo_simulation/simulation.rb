require "thread"

module Rubybayes

  module MonteCarloSimulation
  
    class EngineWrapper
      def initialize(&block)
        @block = block
      end
      def sample
        @block.call
      end
    end
  
    class Simulation
    
      attr_accessor :result
      
      def initialize(&block)
        @results = {}
        #set defaults
        @burn_in  = 0
        @iterations = 1000
        @chains = 1
        @calculate = Proc.new {|x| x}
        @sample = nil
        @engine = nil
        #set user input
        self.instance_eval(&block)
      end
      
      # public functions
      
      def run  
        run_simulation
        get_chains
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
      
      def generate_engine(&block)
        @engine = block
      end
      
      def sample(&block)
        @sample = block
      end
      
      def calculate(&block)
        @calculate = block
      end           
            
      # run simulation  
       
      def run_simulation
        threads = []
        semaphore = Mutex.new
        @chains.times do |nchain|           
          threads << Thread.new do 
            r = run_chain(nchain) 
            semaphore.synchronize { @results[nchain] = r }
          end
        end        
        threads.each {|th| th.join }
      end
       
      def run_chain(nchain)           
        if not (@engine.nil? ^ @sample.nil?)
          raise ArgumentError.new("'sample' or 'generate_sampler' needs to be set for Monte Carlo simulation and not both")
        end
        
        if @engine.nil?
          sampling_engine = EngineWrapper.new(&@sample)
        else
          sampling_engine = @engine.call(nchain)
        end
        
        ret = []
        @burn_in.times {|i| sampling_engine.sample } if @burn_in > 0        
        @iterations.times {|i| ret << @calculate.call( sampling_engine.sample )}
        ret
      end     
          
    end 
    
    # module functions    
    
    def self.merge_chains(results)
      ret = []
      results.each_value {|vals| ret.concat(vals) }
      {0 => ret}    
    end         
     
    def self.run(&block)
      Rubybayes::MonteCarloSimulation::Simulation.new(&block).run
    end
    
    def self.run_and_combine(&block)
      Rubybayes::MonteCarloSimulation::Simulation.new(&block).run.get_merged_chains
    end    
    
  end
  
end