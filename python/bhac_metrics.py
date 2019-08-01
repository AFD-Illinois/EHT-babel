
import numpy as np

def gcov_logks_2d(x1,x2,a, hslope=0, bhac_mks=False):
  """ Returns logks g_munu(x1,x2) for given a.
      Courtesy of (stolen from) George Wong"""

  r = np.exp(x1)
  if bhac_mks:
    th = x2 + 2*hslope/(np.pi**2)*x2*(np.pi - 2*x2)*(np.pi-x2)
  else:
    th = x2

  print(r.shape, th.shape)

  # calculate normal ks
  gcov = np.zeros([x1.shape[0],x1.shape[1],4,4])
  cth = np.cos(th)
  sth = np.sin(th)
  s2 = sth*sth
  rho2 = r*r + a*a*cth*cth
  gcov[:,:,0,0] = (-1. + 2. * r / rho2)
  gcov[:,:,0,1] = (2. * r / rho2)
  gcov[:,:,0,3] = (-2. * a * r * s2 / rho2)
  gcov[:,:,1,0] = gcov[:,:,0,1]
  gcov[:,:,1,1] =  (1. + 2. * r / rho2)
  gcov[:,:,1,3] =  (-a * s2 * (1. + 2. * r / rho2))
  gcov[:,:,2,2] = rho2
  gcov[:,:,3,0] = gcov[:,:,0,3]
  gcov[:,:,3,1] = gcov[:,:,1,3]
  gcov[:,:,3,3] = s2 * (rho2 + a*a * s2 * (1. + 2. * r / rho2))

  if bhac_mks:
    # transform to bhac_mks
    trans = np.zeros([x1.shape[0],x1.shape[1],4,4])
    trans[:,:,0,0] = 1.
    trans[:,:,1,1] = np.exp(x1)
    trans[:,:,2,2] = 1 - 2 * hslope + 12 * hslope * ((x2 / np.pi)**2 - x2/np.pi)
    trans[:,:,3,3] = 1.
  else:
    # transform to logks
    trans = np.zeros([x1.shape[0],x1.shape[1],4,4])
    trans[:,:,0,0] = 1.
    trans[:,:,1,1] = np.exp(x1)
    trans[:,:,2,2] = 1.
    trans[:,:,3,3] = 1.

  #return np.transpose(np.matmul(trans,np.matmul(gcov,trans)),(1,0,2,3))
  return np.matmul(trans,np.matmul(gcov,trans))
