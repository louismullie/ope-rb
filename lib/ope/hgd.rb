module OPE
  
  class HGD
    
    # Random variates from the hypergeometric distribution.
    # Returns the number of white balls drawn when kk balls
    # are drawn at random from an urn containing nn1 white
    # and nn2 black balls.
    
    def self.rhyper(kk, nn1, nn2, coins, precision)

      # "#{coins}, #{precision}"
      
      ix = nil
      
      prng = PRNG.new(coins)
      
      con, deltal = 57.56462733, 0.0078
      deltau, scale = 0.0034, 1.0e25
      
      # Check validity of parameters.
      if [nn1, nn2, kk].include?(Float::INFINITY) ||
         (nn1 < 0 || nn2 < 0 || kk < 0 || kk > (nn1 + nn2))
        raise 'Invalid parameters nn1, nn2 or kk'
      end 
      
      reject = true
      
      if nn1 >= nn2
        n1, n2 = nn2.to_f, nn1.to_f
      else
        n1, n2 = nn1.to_f, nn2.to_f
      end
      
      tn = n1 + n2
      
      if kk + kk >= tn
        k = tn - kk
      else
        k = kk
      end
      
      m = (k+1.0) * (n1+1.0) / (tn+2.0)
      
      minjx = (k - n2 < 0) ? 0 : k - n2
      maxjx = (n1 < k) ? n1 : k
      
      # Degenerate distribution
      if minjx == maxjx
      
        # no, need to untangle TSL
        return maxjx.floor.to_f
      
      # Inverse transformation
      elsif m-minjx < 10
        
        w = nil
        
        if k < n2
          w = afc_op_1(con, n1, n2, k)
        else
          w = afc_op_2(con, n1, n2, k)
        end
        
        count_10 = true
        
        while count_10
          
          count_10 = false
          
          p = w
          ix = minjx
          u = prng.draw * scale
          
          count_20 = true
          
          while count_20
            
            count_20 = false
            
            if u > p
              u = u - p
              p = p * (n1-ix)*(k-ix)
              ix = ix + 1
              p = p / ix / (n2-k+ix)
              if ix > maxjx
                count_10 = true
                break
              end
              count_20 = true
              
            end
            
          end
          
        end
        
      # Hypergeometrics-2 points-exponential tails
      else
        
        s = Math.sqrt( (tn-k) * k * n1 * n2 / (tn - 1.0) / tn / tn)
        
        # Truncation centers cell boundaries at 0.5
        d = (1.5 * s).floor.to_f + 0.5
        xl = m - d + 0.5
        xr = m + d + 0.5
        a = afc_op_3(m, n1, n2, k)
        
        expon = afc_op_4(a, xl, n1, n2, k)
        
        kl = Math.exp(expon)
        
        kr = afc_op_5(a, xr, n1, n2, k)
        
        lamdl = -Math.log(xl *
          (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1))
          
        lamdr = -Math.log(
          (n1 - xr + 1) * (k -xr + 1) / xr / (n2 - k + xr))
        
        p1 = 2 * d
        p2 = p1 + kl / lamdl
        p3 = p2 + kr / lamdr
        
        count_30 = true
        
        while count_30
          
          count_30 = false
          
          u = prng.draw * p3
          v = prng.draw

          # Rectangular region
          if u < p1
            ix = xl + u
          # Left tail region
          elsif u <= p2
            ix = xl + Math.log(v) / lamdl
            if ix < minjx
              count_30 = true
              next 
            end
            v = v * (u - p1) * lamdl
          # Right tail region
          else
            ix = xr - Math.log(v) / lamdr
            if ix > maxjx
              count_30 = true
              next
            end
            v = v * (u - p2) * lamdr
          end
          
          f = nil
          
          if m < 100 || ix <= 50
            
            f = 1
            
            if m < ix
              
              i = m + 1
              while i < ix # <= ?
                f = f * (n1 - i + 1.0) * (k - i + 1.0) / (n2 - k + i) / i
                i += 1
              end
              
            elsif m > ix
              
              i = ix + 1
              while i < m # <= ?
                f = f * i * (n2 - k + i) / (n1 - i) / (k - i) # + 1 ?
                i += 1
              end
              
            end
            
            reject = false if v <= f
            
          else
            
            y = ix; y1 = y + 1.0; ym = y - m
            yn = n1 - y + 1.0; yk = k - y + 1.0
            nk = n2 - k + y1; r = -ym / y1
            s2 = ym / yn; t = ym / yk; e = -ym / nk
            g = yn * yk / (y1 * nk) - 1.0
            dg = 1.0; dg = 1.0 + g if g < 0
            gu = g * (1.0 + g * (-0.5 + g / 3.0))
            gl = gu - 0.25 * (g * g * g * g) / dg
            xm = m + 0.5; xn = n1 - m + 0.5
            xk = k - m + 0.5; nm = n2 - k + xm
            
            ub = y * gu - m * gl + deltau + xm *
                 r * (1 + r * (-0.5 + r / 3)) +
                 xn * s2 * (1.0 + s2 * (-0.5 + s2 / 3.0)) +
                 xk * t * (1.0 + t * (-0.5 + t / 3.0)) +
                 nm * e * (1.0 + e * (-0.5 + e / 3.0))
            
            alv = Math.log(v)
            
            if alv > ub
              reject = true
            else
              
              dr = xm * (r * r * r * r)
              dr /= (1.0 + r) if r < 0
              ds = xn * (s2 * s2 * s2 * s2)
              ds /= (1.0 + s2) if s2 < 0
              dt = xk * (t * t * t * t)
              dt /= (1.0 + t) if t < 0
              de = nm * (e * e * e * e)
              de /= (1.0 + e) if e < 0
              
              
              cand = ub - 0.25 * (dr + ds + dt + de) +
              (y + m) * (gl - gu) - deltal
            
              if alv < cand
                reject = false
              else  
              
            
              cand = afc_op_6(a, ix, n1, n2, k)
              
              if alv <= cand
                reject = false
              else
                reject = true
              end
            
            end
          
          end
          
        end
        
        count_30 = true if reject
        
      end
        
      end
      
      if kk + kk >= tn
        if nn1 > nn2
          ix = kk - nn2 + ix
        else
          ix = nn1 - ix
        end
      else
        ix = kk - ix if nn1 > nn2
      end
      
      jx = ix
      
      return jx.floor.to_f
      
    end
    
    def self.afc_op_1(con, n1, n2, k)
      Math.exp(con + afc(n2) + afc(n1+n2-k) - afc(n2-k) - afc(n1+n2))
    end
    
    def self.afc_op_2(con, n1, n2, k)
      Math.exp(con + afc(n1) + afc(k) + afc(k-n2) - afc(n1+n2))
    end
    
    def self.afc_op_3(m, n1, n2, k)
      afc_native_op_3(m, n1, n2, k).to_f
    end
    
    def self.afc_op_4(a, xl, n1, n2, k)
      a - afc(xl) - afc(n1 - xl) - afc(k - xl) - afc(n2 - k + xl)
    end
    
    def self.afc_op_5(a, xr, n1, n2, k)
      afc_native_op_5(a, xr, n1, n2, k)
    end
    
    def self.afc_op_6(a, ix, n1, n2, k)
      a - afc(ix) - afc(n1 - ix) - afc(k - ix) - afc(n2 - k + ix)
    end
    
    def self.afc(i)
      
      afc_native(i).to_f
      
    end
    
    
  
  end
  
end