# Following code is based on:
#
# Minimization- Minimization algorithms on pure Ruby
# Copyright (C) 2010 Claudio Bustos
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

module Rubybayes

  module Minimization
    FailedIteration=Class.new(Exception)
    # Base class for unidimensional minimizers
    class Unidimensional
      # Default value for error on f(x)
      EPSILON=1e-6
      # Default number of maximum iterations
      MAX_ITERATIONS=100
      # Minimum value for x
      attr_reader :x_minimum
      # Minimum value for f(x)
      attr_reader :f_minimum
      # Log of iterations. Should be an array
      attr_reader :log
      # Name of fields of log
      attr_reader :log_header
      # Absolute error on x
      attr_accessor :epsilon
      # Expected value. Fast minimum finding if set
      attr_reader :expected
      # Numbers of iterations
      attr_reader :iterations
      # Create a new minimizer
      def initialize(lower, upper, proc)
        raise "first argument  should be lower than second" if lower>=upper
        @lower=lower
        @upper=upper
        @proc=proc
        golden = 0.3819660;
        @expected = @lower + golden * (@upper - @lower);
        @max_iteration=MAX_ITERATIONS
        @epsilon=EPSILON
        @iterations=0
        @log=[]
        @log_header=%w{I xl xh f(xl) f(xh) dx df(x)}
      end
      # Set expected value
      def expected=(v)
        @expected=v
      end
      def log_summary
        @log.join("\n")
      end
      # Convenience method to minimize
      # == Parameters:
      # * <tt>lower</tt>: Lower possible value
      # * <tt>upper</tt>: Higher possible value
      # * <tt>expected</tt>: Optional expected value. Faster the search is near correct value.
      # * <tt>&block</tt>: Block with function to minimize
      # == Usage:
      #   minimizer=Minimization::GoldenSection.minimize(-1000, 1000) {|x|
      #             x**2 }
      # 
      def self.minimize(lower,upper,expected=nil,&block)
        minimizer=new(lower,upper,block)
        minimizer.expected=expected unless expected.nil?
        raise FailedIteration unless minimizer.iterate
        minimizer
      end
      # Iterate to find the minimum
      def iterate
        raise "You should implement this"
      end
      def f(x)
        @proc.call(x)
      end
    end
    
    class Brent < Unidimensional
      GSL_SQRT_DBL_EPSILON=1.4901161193847656e-08
      def initialize(lower,upper, proc)
        super

        @do_bracketing=true

        # Init

        golden = 0.3819660;      #golden = (3 - sqrt(5))/2

        v = @lower + golden * (@upper - @lower);
        w = v;

        @x_minimum = v ;
        @f_minimum = f(v) ;
        @x_lower=@lower
        @x_upper=@upper
        @f_lower = f(@lower) ;
        @f_upper = f(@lower) ;

        @v = v;
        @w = w;

        @d = 0;
        @e = 0;
        @f_v=f(v)
        @f_w=@f_v
      end

      def expected=(v)
        @x_minimum=v
        @f_minimum=f(v)
        @do_bracketing=false
      end
      
      def bracketing
        eval_max=10
        f_left = @f_lower;
        f_right = @f_upper;
        x_left = @x_lower;
        x_right= @x_upper;
        golden = 0.3819660;      # golden = (3 - sqrt(5))/2 */
        nb_eval=0

        if (f_right >= f_left)
          x_center = (x_right - x_left) * golden + x_left;
          nb_eval+=1;
          f_center=f(x_center)
        else
          x_center = x_right ;
          f_center = f_right ;
          x_right = (x_center - x_left).quo(golden) + x_left;
          nb_eval+=1;
          f_right=f(x_right);
        end


        begin
          @log << ["B#{nb_eval}", x_left, x_right, f_left, f_right, (x_left-x_right).abs, (f_left-f_right).abs]
          if (f_center < f_left )
            if (f_center < f_right)
              @x_lower = x_left;
              @x_upper = x_right;
              @x_minimum = x_center;
              @f_lower = f_left;
              @f_upper = f_right;
              @f_minimum = f_center;
              return true;
            elsif (f_center > f_right)
              x_left = x_center;
              f_left = f_center;
              x_center = x_right;
              f_center = f_right;
              x_right = (x_center - x_left).quo(golden) + x_left;
              nb_eval+=1;
              f_right=f(x_right);
            else # f_center == f_right */
              x_right = x_center;
              f_right = f_center;
              x_center = (x_right - x_left).quo(golden) + x_left;
              nb_eval+=1;
              f_center=f(x_center);
            end
          else # f_center >= f_left */
            x_right = x_center;
            f_right = f_center;
            x_center = (x_right - x_left) * golden + x_left;
            nb_eval+=1;
            f_center=f(x_center);
          end
        end while ((nb_eval < eval_max) and
        ((x_right - x_left) > GSL_SQRT_DBL_EPSILON * ( (x_right + x_left) * 0.5 ) + GSL_SQRT_DBL_EPSILON))
        @x_lower = x_left;
        @x_upper = x_right;
        @x_minimum = x_center;
        @f_lower = f_left;
        @f_upper = f_right;
        @f_minimum = f_center;
        return false;

      end
      # Start the minimization process
      # If you want to control manually the process, use brent_iterate
      def iterate
        k=0
        bracketing if @do_bracketing
        while k<@max_iteration and (@x_lower-@x_upper).abs>@epsilon
          k+=1
          result=brent_iterate
          raise FailedIteration,"Error on iteration" if !result
          begin 
            @log << [k, @x_lower, @x_upper, @f_lower, @f_upper, (@x_lower-@x_upper).abs, (@f_lower-@f_upper).abs]
          rescue =>@e
            @log << [k, @e.to_s,nil,nil,nil,nil,nil]
          end
        end
        @iterations=k
        return true
      end
      # Generate one iteration.
      def brent_iterate
        x_left = @x_lower;
        x_right = @x_upper;

        z = @x_minimum;
        d = @e;
        e = @d;
        v = @v;
        w = @w;
        f_v = @f_v;
        f_w = @f_w;
        f_z = @f_minimum;

        golden = 0.3819660;      # golden = (3 - sqrt(5))/2 */

        w_lower = (z - x_left)
        w_upper = (x_right - z)

        tolerance =  GSL_SQRT_DBL_EPSILON * z.abs

        midpoint = 0.5 * (x_left + x_right)
        _p,q,r=0,0,0
        if (e.abs > tolerance)

          # fit parabola */

          r = (z - w) * (f_z - f_v);
          q = (z - v) * (f_z - f_w);
          _p = (z - v) * q - (z - w) * r;
          q = 2 * (q - r);

          if (q > 0)
            _p = -_p
          else
            q = -q;
          end
          r = e;
          e = d;
        end

        if (_p.abs < (0.5 * q * r).abs and _p < q * w_lower and _p < q * w_upper)
          t2 = 2 * tolerance ;

          d = _p.quo(q);
          u = z + d;

          if ((u - x_left) < t2 or (x_right - u) < t2)
            d = (z < midpoint) ? tolerance : -tolerance ;
          end
        else

          e = (z < midpoint) ? x_right - z : -(z - x_left) ;
          d = golden * e;
        end

        if ( d.abs >= tolerance)
          u = z + d;
        else
          u = z + ((d > 0) ? tolerance : -tolerance) ;
        end

        @e = e;
        @d = d;

        f_u=f(u)

        if (f_u <= f_z)
          if (u < z)
            @x_upper = z;
            @f_upper = f_z;
          else
            @x_lower = z;
            @f_lower = f_z;
          end
          @v = w;
          @f_v = f_w;
          @w = z;
          @f_w = f_z;
          @x_minimum = u;
          @f_minimum = f_u;
          return true;
        else
          if (u < z)
            @x_lower = u;
            @f_lower = f_u;
            return true;
          else
            @x_upper = u;
            @f_upper = f_u;
            return true;
          end

          if (f_u <= f_w or w == z)
            @v = w;
            @f_v = f_w;
            @w = u;
            @f_w = f_u;
            return true;
          elsif f_u <= f_v or v == z or v == w
            @v = u;
            @f_v = f_u;
            return true;
          end

        end
        return false

      end
      
    end
    
  end
  
  module Minimization

    # class which holds the point,value pair
    class PointValuePair      
      attr_accessor :value
      attr_reader   :point

      # == Parameters:
      # * <tt>point</tt>: Coordinates of the point
      # * <tt>value</tt>: Function value at the point
      #
      def initialize(point, value)
        @point = point.clone
        @value  = value
      end

      # returns a copy of the point
      def get_point_clone
        return @point.clone
      end
    end

    class DirectSearchMinimizer

      EPSILON_DEFAULT         = 1e-6
      MAX_ITERATIONS_DEFAULT  = 1000000

      attr_reader :x_minimum
      attr_reader :f_minimum
      attr_reader :epsilon

      def initialize(f, start_point, iterate_simplex_ref)
        @epsilon             = EPSILON_DEFAULT
        # Default number of maximum iterations
        @max_iterations      = MAX_ITERATIONS_DEFAULT
        # proc which iterates the simplex
        @iterate_simplex_ref = iterate_simplex_ref
        @relative_threshold  = 100 * @epsilon
        @absolute_threshold  = @epsilon
        @x_minimum           = nil
        @f_minimum           = nil
        @f = f

        # create and initializ start configurations
        if @start_configuration == nil
          # sets the start configuration point as unit
          self.start_configuration = Array.new(start_point.length) { 1.0 }
        end

        @iterations  = 0
        @evaluations = 0
        # create the simplex for the first time
        build_simplex(start_point)
        evaluate_simplex
      end

      def f(x)
        return @f.call(x)
      end

      def iterate_simplex
        return iterate_simplex_ref.call
      end

      # increment iteration counter by 1
      def increment_iterations_counter
        @iterations += 1
        raise "iteration limit reached" if @iterations > @max_iterations
      end

      # compares 2 PointValuePair points
      def compare(v1, v2)
        if v1.value == v2.value
          return 0
        elsif v1.value > v2.value
          return 1
        else
          return -1
        end
      end

      # checks whether the function is converging
      def converging?
        # check the convergence in a given direction comparing the previous and current values
        def point_converged?(previous, current)
          pre        = previous.value
          curr       = current.value
          diff       = (pre - curr).abs
          size       = [pre.abs, curr.abs].max
          return !((diff <= (size * @relative_threshold)) and (diff <= @absolute_threshold))
        end

        # returns true if converging is possible atleast in one direction
        if @iterations > 0
          # given direction is converged
          converged = true
          0.upto(@simplex.length - 1) do |i|
            converged &= !point_converged?(@previous[i], @simplex[i])
          end
          return !converged
        end

        # if no iterations were done, convergence undefined
        return true
      end

      # only the relative position of the n vertices with respect
      # to the first one are stored
      def start_configuration=(steps)
        n = steps.length
        @start_configuration = Array.new(n) { Array.new(n, 0) }
        0.upto(n - 1) do |i|
          vertex_i = @start_configuration[i]
          0.upto(i) do |j|
            raise "equals vertices #{j} and #{j+1} in simplex configuration" if steps[j] == 0.0
            0.upto(j) do |k|
              vertex_i[k] = steps[k]
            end
          end
        end
      end

      # Build an initial simplex
      # == Parameters:
      # * <tt>start_point</tt>: starting point of the minimization search
      #
      def build_simplex(start_point)
        n = start_point.length
        raise "dimension mismatch" if n != @start_configuration.length
        # set first vertex
        @simplex = Array.new(n+1)
        @simplex[0] = PointValuePair.new(start_point, Float::NAN)

        # set remaining vertices
        0.upto(n - 1) do |i|
          conf_i   = @start_configuration[i]
          vertex_i = Array.new(n)
          0.upto(n - 1) do |k|
            vertex_i[k] = start_point[k] + conf_i[k]
          end
          @simplex[i + 1] = PointValuePair.new(vertex_i, Float::NAN)
        end
      end

      # Evaluate all the non-evaluated points of the simplex
      def evaluate_simplex
        # evaluate the objective function at all non-evaluated simplex points
        0.upto(@simplex.length - 1) do |i|
          vertex = @simplex[i]
          point  = vertex.point
          if vertex.value.nan?
            @simplex[i] = PointValuePair.new(point, f(point))
          end
        end
        # sort the simplex from best to worst
        @simplex.sort!{ |x1, x2| x1.value <=> x2.value }
      end

      # Replace the worst point of the simplex by a new point
      # == Parameters:
      # * <tt>point_value_pair</tt>: point to insert
      #
      def replace_worst_point(point_value_pair)
        n = @simplex.length - 1
        0.upto(n - 1) do |i|
          if (compare(@simplex[i], point_value_pair) > 0)
            point_value_pair, @simplex[i] = @simplex[i], point_value_pair
          end
        end
        @simplex[n] = point_value_pair
      end

      # Convenience method to minimize
      # == Parameters:
      # * <tt>start_point</tt>: Starting points
      # * <tt>f</tt>: Function to minimize
      # == Usage:
      #   minimizer=Minimization::NelderMead.minimize(proc{|x| (x[0] - 1) ** 2 + (x[1] - 5) ** 2}, [0, 0])
      #
      def self.minimize(f, start_point)
        min=Minimization::NelderMead.new(f, start_point)
        while min.converging?
          min.iterate
        end
        return min
      end

      # Iterate the simplex one step. Use this when iteration needs to be done manually
      # == Usage:
      #   minimizer=Minimization::NelderMead.new(proc{|x| (x[0] - 1) ** 2 + (x[1] - 5) ** 2}, [0, 0])
      #   while minimizer.converging?
      #     minimizer.Iterate
      #   end
      #   minimizer.x_minimum
      #   minimizer.f_minimum
      #
      def iterate
        # set previous simplex as the current simplex
        @previous = Array.new(@simplex.length)
        0.upto(@simplex.length - 1) do |i|
          point = @simplex[i].point                                # clone require?
          @previous[i] = PointValuePair.new(point, f(point))
        end
        # iterate simplex
        iterate_simplex
        # set results
        @x_minimum = @simplex[0].point
        @f_minimum = @simplex[0].value
      end
    end

    # = Nelder Mead Minimizer.
    # A multidimensional minimization methods.
    # == Usage.
    #  require 'minimization'
    #  min=Minimization::NelderMead.new(proc {|x| (x[0] - 2)**2 + (x[1] - 5)**2}, [1, 2])
    #  while min.converging?
    #    min.iterate
    #  end
    #  min.x_minimum
    #  min.f_minimum
    #
    class NelderMead < DirectSearchMinimizer
      def initialize(f, start_point)
        # Reflection coefficient
        @rho   = 1.0
        # Expansion coefficient
        @khi   = 2.0
        # Contraction coefficient
        @gamma = 0.5
        # Shrinkage coefficient
        @sigma = 0.5
        super(f, start_point, proc{iterate_simplex})
      end
      
      def iterate_simplex
        increment_iterations_counter
        n = @simplex.length - 1
        # the simplex has n+1 point if dimension is n
        best       = @simplex[0]
        secondBest = @simplex[n - 1]
        worst      = @simplex[n]
        x_worst    = worst.point
        centroid = Array.new(n, 0)
        # compute the centroid of the best vertices
        # (dismissing the worst point at index n)
        0.upto(n - 1) do |i|
          x = @simplex[i].point
          0.upto(n - 1) do |j|
            centroid[j] += x[j]
          end
        end
        scaling = 1.0 / n
        0.upto(n - 1) do |j|
          centroid[j] *= scaling
        end
        xr = Array.new(n)
        # compute the reflection point
        0.upto(n - 1) do |j|
          xr[j] = centroid[j] + @rho * (centroid[j] - x_worst[j])
        end
        reflected = PointValuePair.new(xr, f(xr))
        if ((compare(best, reflected) <= 0) && (compare(reflected, secondBest) < 0))
          # accept the reflected point
          replace_worst_point(reflected)
        elsif (compare(reflected, best) < 0)
          xe = Array.new(n)
          # compute the expansion point
          0.upto(n - 1) do |j|
            xe[j] = centroid[j] + @khi * (xr[j] - centroid[j])
          end
          expanded = PointValuePair.new(xe, f(xe))
          if (compare(expanded, reflected) < 0)
            # accept the expansion point
            replace_worst_point(expanded)
          else
            # accept the reflected point
            replace_worst_point(reflected)
          end
        else
          if (compare(reflected, worst) < 0)
            xc = Array.new(n)
            # perform an outside contraction
            0.upto(n - 1) do |j|
              xc[j] = centroid[j] + @gamma * (xr[j] - centroid[j])
            end
            out_contracted = PointValuePair.new(xc, f(xc))
            if (compare(out_contracted, reflected) <= 0)
              # accept the contraction point
              replace_worst_point(out_contracted)
              return
            end
          else
            xc = Array.new(n)
            # perform an inside contraction
            0.upto(n - 1) do |j|
              xc[j] = centroid[j] - @gamma * (centroid[j] - x_worst[j])
            end
            in_contracted = PointValuePair.new(xc, f(xc))

            if (compare(in_contracted, worst) < 0)
              # accept the contraction point
              replace_worst_point(in_contracted)
              return
            end
          end
          # perform a shrink
          x_smallest = @simplex[0].point
          0.upto(@simplex.length - 1) do |i|
            x = @simplex[i].get_point_clone
            0.upto(n - 1) do |j|
              x[j] = x_smallest[j] + @sigma * (x[j] - x_smallest[j])
            end
            @simplex[i] = PointValuePair.new(x, Float::NAN)
          end
          evaluate_simplex
        end
      end
    end
  end  

end