module Rubybayes

  module MonteCarloSimulation
  
    class Measurement
      attr_reader :data
      def initialize(vector = [])
        @data = vector
        @sort = nil
      end

      def method_missing(name, *args, &block)
        @data.send(name, *args, &block)
      end
	
      def print_stats
        puts "Sum: #{sum.round(2)}, Mean: #{mean.round(2)}, Min: #{min.round(2)}, Max: #{max.round(2)}"
      end
      
      def sort
        @sort ||= @data.sort
      end
    
      def sum
        @data.inject(0.0) {|sum,x| sum + x }
      end
	
      def mean
        sum / size
      end
      
      def variance #sample variance
        @data.inject(0.0) {|sum,x| sum + (x - mean)**2} / (size - 1)
      end
      
      def sd
        Math.sqrt(variance)
      end
      
      def quantile(p)
        idx = (p.to_f * (size - 1)).round(6)
        upper = idx.ceil
        lower = idx.floor        
        if upper == lower
          sort[lower]
        else
          frac = idx - lower  
          frac * sort[upper] + (1.0 - frac) * sort[lower]
        end
      end
      
      def cov(b)
        raise "cov error - dimensions of measurements don't match" if size != b.size
        sum = 0.0
        size.times do |i|
          sum += (@data[i] - mean) * (b[i] - b.mean)
        end
        sum / (size - 1)
      end
      
      def cor(b)        
        cov(b) / (sd * b.sd)
      end
      
      def acf(lag = 0)
        raise ArgumentError.new("acf error - lag probably too big") if lag >= size
        x = @data.take(size-lag)
        y = @data.rotate(lag).take(size-lag)
        raise "acf error - problem calculating lagged vectors" if x.size != y.size
        sum = 0.0
        (size-lag).times do |i|
          sum += (x[i] - mean) * (y[i] - mean)
        end
        sum / ((size - 1) * variance)
      end	
     
      def cumulative_error
        squared = @data.each_with_index.map {|h, i| (h - cumulative_mean[i])**2 }
        cum = create_cumulative(squared).map{|x| Math.sqrt(x)}
        div(cum,cumulative_size)
      end
	
      def cumulative_size
        (1..@data.size).to_a		
      end
	
      def cumulative_mean
        div(cumulative_sum,cumulative_size)
      end
	
      def cumulative_sum
        create_cumulative(@data)
      end
	
      private
            
      def div(a,b)
        a.zip(b).map{|x, y| x / y}
      end
	
      def create_cumulative(vec)
        sum = 0.0
        vec.map{|x| sum += x}
      end
      
    end

    def self.extract_measurements(result)      
      # results as hash {1 => [[A,B,C],[A,B,C],..], 2 => [[A,B,C],[A,B,C],..], ..}
      # transposes into hash {1 => [[A,A,..],[B,B,..],[C,C,..]], 2 => ..}
      # and returns hash  {1 => [Measurement.new([A,A,..]), Measurement.new([B,B,..]), ...], 2=> ..}
      ret = {}
      result.each_pair do |key, value| 
        next if value.size < 1
        if value[0].kind_of?(Array)
          ret[key] = value.transpose.collect {|x| Measurement.new(x) }
        else
          ret[key] = Measurement.new(value)
        end
      end      
      ret
    end
    
  end
    
end