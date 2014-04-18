module OPE
  
  class HGD
    
    AFCTable = [
      '0.0', # ln(0!) = ln(1)
      '0.0', # ln(1!) = ln(1)
      '0.69314718055994530941723212145817', # ln(2!)
      '1.79175946922805500081247735838070', # ln(3!)
      '3.17805383034794561964694160129705', # ln(4!)
      '4.78749174278204599424770093452324', # ln(5!)
      '6.57925121201010099506017829290394', # ln(6!)
      '8.52516136106541430016553103634712', # ln(7!)
      '10.60460290274525022841722740072165' # ln(8!)
    ].map { |x| BigDecimal.new(x) }
    
    # Random variates from the hypergeometric distribution.
    # Returns the number of white balls drawn when kk balls
    # are drawn at random from an urn containing nn1 white
    # and nn2 black balls.
    
    def self.rhyper(kk, nn1, nn2, coins, precision)

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
          w = Math.exp(con + afc(n2) + afc(n1+n2-k)-afc(n2-k)-afc(n1+n2))
        else
          w = Math.exp(con + afc(n1) + afc(k) + afc(k-n2) -afc(n1+n2))
        end
        
        catch :l10 do
          
          p = w
          ix = minjx
          u = prng.draw * scale
          
          catch :l20 do
            
            if u > p
              u = u - p
              p = p * (n1-ix)*(k-ix)
              ix = ix + 1
              p = p / ix / (n2-k+ix)
              throw :l10 if ix > maxjx
              throw :l20
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
        a = afc(m) + afc(n1-m) +
          afc(k-m) + afc(n2-k+m) 
        
        expon = a - afc(xl) -
          afc(n1 - xl) - afc(k - xl) -
          afc(n2 - k + xl)
        
        kl = Math.exp(expon)
        
        kr = Math.exp(a - afc(xr-1) -
          afc(n1-xr+1) - afc(k-xr+1) -
          afc(n2-k+xr-1))
        
        lamdl = -Math.log(xl *
          (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1))
          
        lamdr = -Math.log(
          (n1 - xr + 1) * (k -xr + 1) / xr / (n2 - k + xr))
        
        p1 = 2 * d
        p2 = p1 + kl / lamdl
        p3 = p2 + kr / lamdr
        
        count_30 = 0
        
        catch :count do
          
          count_30 += 1

          u = prng.draw * p3
          v = prng.draw

          # Rectangular region
          if u < p1
            ix = xl + u
          # Left tail region
          elsif u <= p2
            ix = xl + Math.log(v) / lamdl
            throw :count if ix < minjx
            v = v * (u - p1) * lamdl
          # Right tail region
          else
            ix = xr - Math.log(v) / lamdr
            throw :count if ix > maxjx
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
              
            
              cand = a - afc(ix) - afc(n1 - ix) - afc(k - ix) - afc(n2 - k + ix)
              
              if alv <= cand
                reject = false
              else
                reject = true
              end
            
            end
          
          end
          
        end
        
        throw :count if reject
        
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
    
    # Calculates logarithm of i factorial: ln(i!)
    # If i < 9, uses table lookup. Otherwise, uses
    # Stirling's approximation.
    def self.afc(i)
      
      raise 'i should not be < 0' if i < 0
      
      return AFCTable[i] if i <= 8
      
      frac_12, frac_360 = 1.0 / 12.0, 1.0 / 360.0
      frac_pi = 0.5 * Math.log(2 * Math::PI)
      
      (i + 0.5) * Math.log(i) - i + frac_12 /
      i - frac_360 / i / i / i + frac_pi
      
    end
  
  end
  
end