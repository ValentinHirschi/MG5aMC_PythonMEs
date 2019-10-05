##########################################################################################
#
# Copyright (c) 2017 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
##########################################################################################

import logging
import math
import copy
import numpy as np

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')


class InvalidOperation(Exception):
    pass

def almost_equal(x, y, rel_tol=0, abs_tol=0):
    """Check if two objects are equal within certain relative and absolute tolerances.
    The operations abs(x + y) and abs(x - y) need to be well-defined
    for this function to work.

    :param x: first object to be compared
    :param y: second object to be compared
    :param rel_tol: relative tolerance for the comparison
    :param abs_tol: absolute tolerance for the comparison
    :return: True if the elements are equal within tolerance, False otherwise
    :rtype: bool
    """

    diffxy = abs(x - y)
    if diffxy <= abs_tol: return True
    sumxy = abs(x + y)
    # Rough check that the ratio is smaller than 1 to avoid division by zero
    if sumxy < diffxy: return False
    return diffxy / sumxy <= rel_tol

#=========================================================================================
# Vector
#=========================================================================================

class Vector(np.ndarray):

    def __new__(cls, *args, **opts):

        if args and isinstance(args[0], Vector):
            foo = args[0].get_copy()
        else:
            foo = np.asanyarray(*args, **opts).view(cls)
        return foo

    def __array_finalize__(self, obj):

        try: self.eps = np.finfo(self.dtype).eps ** 0.5
        except: self.eps = 0
        return

    def huge(self):

        if np.issubdtype(self.dtype, np.inexact):
            return np.finfo(self.dtype).max
        elif np.issubdtype(self.dtype, np.integer):
            return np.iinfo(self.dtype).max
        else:
            raise ValueError

    def almost_zero(self, x):

        return x < self.eps

    def almost_equal(self, x, y=None):
        """Check if two numbers are equal within the square root
        of the typical numerical accuracy of the underlying array data type.
        """

        if y is None: y = self
        return almost_equal(x, y, rel_tol=self.eps)

    def square(self):

        return self.dot(self)

    def __abs__(self):

        foo = self.view(np.ndarray)
        return np.dot(foo, foo) ** 0.5

    def __eq__(self, other):

        return almost_equal(self, other, rel_tol=self.eps+other.eps)

    def __hash__(self):

        return tuple(x for x in self).__hash__()

    def get_copy(self):

        # The vector instantiated by get_copy() should be modified
        # without changing the previous instance, irrespectively of the
        # (presumably few) layers that compose entries of the vector
        # return copy.deepcopy(self)
        return copy.copy(self)

    def normalize(self):

        self.__idiv__(abs(self))
        return self

    def project_onto(self, v):

        return (self.dot(v) / v.square()) * v

    def component_orthogonal_to(self, v):

        return self - self.project_onto(v)

    @classmethod
    def cos(cls, v, w):
        """Cosine of the angle between two vectors."""

        assert v.square() > 0
        assert w.square() > 0
        return v.dot(w)/(abs(v)*abs(w))

    # Specific to 3D vectors
    def cross(self, v):

        assert len(self) == 3
        assert len(v) == 3
        return Vector([
            self[1] * v[2] - self[2] * v[1],
            self[2] * v[0] - self[0] * v[2],
            self[0] * v[1] - self[1] * v[0]
        ])

#=========================================================================================
# LorentzVector
#=========================================================================================

class LorentzVector(Vector):

    def __new__(cls, *args, **opts):

        if len(args) == 0:
            return super(LorentzVector, cls).__new__(cls, [0., 0., 0., 0.], **opts)
        return super(LorentzVector, cls).__new__(cls, *args, **opts)

    def space(self):
        """Return the spatial part of this LorentzVector."""

        return self[1:].view(type=Vector)

    def dot(self, v, out=None):
        """Compute the Lorentz scalar product."""
        ## The implementation below allows for a check but it should be done upstream and
        ## significantly slows down the code here.
        # pos = self[0]*v[0]
        # neg = self.space().dot(v.space())
        # if pos+neg != 0 and abs(2*(pos-neg)/(pos+neg)) < 100.*self.eps(): return 0
        # return pos - neg
        return self[0]*v[0] - self[1]*v[1] - self[2]*v[2] - self[3]*v[3]

    def square_almost_zero(self):
        """Check if the square of this LorentzVector is zero within numerical accuracy."""

        return self.almost_zero(self.square() / np.dot(self, self))

    def rho2(self):
        """Compute the radius squared."""

        return self.space().square()

    def rho(self):
        """Compute the radius."""

        return abs(self.space())

    def space_direction(self):
        """Compute the corresponding unit vector in ordinary space."""

        return self.space()/self.rho()

    def set_square(self, square, negative=False):
        """Change the time component of this LorentzVector
        in such a way that self.square() = square.
        If negative is True, set the time component to be negative,
        else assume it is positive.
        """

        # Note: square = self[0]**2 - self.rho2(),
        # so if (self.rho2() + square) is negative, self[0] is imaginary.
        # Letting math.sqrt fail if data is not complex on purpose in this case.
        self[0] = (self.rho2() + square) ** 0.5
        if negative: self[0] *= -1
        return self

    def rotoboost(self, p, q):
        """Apply the Lorentz transformation that sends p in q to this vector."""

        # NOTE: when applying the same Lorentz transformation to many vectors,
        #       this function goes many times through the same checks.

        # Compute squares
        p2 = p.square()
        q2 = q.square()
        # Check if both Lorentz squares are small compared to the euclidean squares,
        # in which case the alternative formula should be used
        if p.square_almost_zero() and q.square_almost_zero():
            # Use alternative formula
            if self.almost_equal(p):
                for i in range(len(self)):
                    self[i] = q[i]
            else:
                logger.critical("Error in vectors.rotoboost: missing formula")
                logger.critical("Boosting %s (%.9e)" % (str(self), self.square()))
                logger.critical("p = %s (%.9e)" % (str(p), p2))
                logger.critical("q = %s (%.9e)" % (str(q), q2))
                logger.critical("Eq. (4.14) of arXiv:0706.0017v2, p. 26 not implemented")
                raise NotImplementedError
            return self
        else:
            # Check that the two invariants are close,
            # else the transformation is invalid
            if not almost_equal(p2, q2, rel_tol=p.eps+q.eps):
                logger.critical("Error in vectors.rotoboost: nonzero, unequal squares")
                logger.critical("p = %s (%.9e)" % (str(p), p2))
                logger.critical("q = %s (%.9e)" % (str(q), q2))
                print "Error in vectors.rotoboost: nonzero, unequal squares"
                print "p = %s (%.9e)" % (str(p), p2)
                print "q = %s (%.9e)" % (str(q), q2)
                raise InvalidOperation
            # Compute scalar products
            pq = p + q
            pq2 = pq.square()
            p_s = self.dot(p)
            pq_s = self.dot(pq)
            # Assemble vector
            self.__iadd__(2 * ((p_s/q2) * q - (pq_s/pq2) * pq))
            return self

    def pt(self, axis=3):
        """Compute transverse momentum."""

        return math.sqrt(
            sum(self[i]**2 for i in range(1, len(self)) if i != axis) )

    def pseudoRap(self):
        """Compute pseudorapidity."""

        pt = self.pt()
        if pt < self.eps and abs(self[3]) < self.eps:
            return self.huge()*(self[3]/abs(self[3]))
        th = math.atan2(pt, self[3])
        return -math.log(math.tan(th/2.))

    def rap(self):
        """Compute rapidity in the lab frame. (needs checking)"""

        if self.pt() < self.eps and abs(self[3]) < self.eps:
            return self.huge()*(self[3]/abs(self[3]))

        return .5*math.log((self[0]+self[3])/(self[0]-self[3]))

    def getdelphi(self, p2):
        """Compute the phi-angle separation with p2."""

        pt1 = self.pt()
        pt2 = p2.pt()
        if pt1 == 0. or pt2 == 0.:
            return self.huge()
        tmp = self[1]*p2[1] + self[2]*p2[2]
        tmp /= (pt1*pt2)
        if abs(tmp) > (1.0+self.eps):
            logger.critical("Cosine larger than 1. in phase-space cuts.")
            raise ValueError
        if abs(tmp) > 1.0:
            return math.acos(tmp/abs(tmp))
        return math.acos(tmp)

    def deltaR(self, p2):
        """Compute the deltaR separation with momentum p2."""

        delta_eta = self.pseudoRap() - p2.pseudoRap()
        delta_phi = self.getdelphi(p2)
        return math.sqrt(delta_eta**2 + delta_phi**2)

    def boostVector(self):

        if self == LorentzVector():
            return Vector([0.] * 3)
        if self[0] <= 0. or self.square() < 0.:
            logger.critical("Attempting to compute a boost vector from")
            logger.critical("%s (%.9e)" % (str(self), self.square()))
            raise InvalidOperation
        return self.space()/self[0]

    def cosTheta(self):

        ptot = self.rho()
        assert (ptot > 0.)
        return self[3] / ptot

    @classmethod
    def cos(cls, v, w):
        """Cosine of the angle between the space part of two vectors."""

        return Vector.cos(v.space(), w.space())

    def phi(self):

        return math.atan2(self[2], self[1])
    
    def boost(self, boost_vector, gamma=-1.):
        """Transport self into the rest frame of the boost_vector in argument.
        This means that the following command, for any vector p=(E, px, py, pz)
            p.boost(-p.boostVector())
        transforms p to (M,0,0,0).
        """

        b2 = boost_vector.square()
        if gamma < 0.:
            gamma = 1.0 / math.sqrt(1.0 - b2)

        bp = self.space().dot(boost_vector)
        gamma2 = (gamma-1.0) / b2 if b2 > 0 else 0.
        factor = gamma2*bp + gamma*self[0]
        self_space = self.space()
        self_space += factor*boost_vector
        self[0] = gamma*(self[0] + bp)
        return self

    @classmethod
    def boost_vector_from_to(cls, p, q):
        """Determine the boost vector for a pure boost that sends p into q.
        For details, see appendix A.2.2 of Simone Lionetti's PhD thesis.

        :param LorentzVector p: Starting Lorentz vector to define the boost.
        :param LorentzVector q: Target Lorentz vector to define the boost.
        :return: Velocity vector for a boost that sends p into q.
        :rtype: Vector
        """

        eps = p.eps+q.eps
        p_abs = abs(p)
        q_abs = abs(q)
        assert almost_equal(p.square(), q.square(), rel_tol=eps) or \
               (p.square_almost_zero() and q.square_almost_zero())
        p_vec = p.space()
        q_vec = q.space()
        if almost_equal(p_vec, q_vec, rel_tol=eps):
            return Vector([0 for _ in p_vec])
        n_vec = (q_vec - p_vec).normalize()
        na = LorentzVector([1, ] + list(+n_vec))
        nb = LorentzVector([1, ] + list(-n_vec))
        assert na.square_almost_zero()
        assert nb.square_almost_zero()
        assert almost_equal(na.dot(nb), 2, rel_tol=eps)
        p_plus  =  p.dot(nb)
        p_minus =  p.dot(na)
        q_plus  =  q.dot(nb)
        q_minus =  q.dot(na)
        if p_minus/p_abs < eps and q_minus/q_abs < eps:
            if p_plus/p_abs < eps and q_plus/q_abs < eps:
                exppy = 1
            else:
                exppy = q_plus / p_plus
        else:
            if p_plus/p_abs < eps and q_plus/q_abs < eps:
                exppy = p_minus / q_minus
            else:
                exppy = ((q_plus*p_minus) / (q_minus*p_plus)) ** 0.5
        expmy = 1. / exppy
        return abs((exppy - expmy) / (exppy + expmy)) * n_vec

    def boost_from_to(self, p, q):
        """Apply a pure boost that sends p into q to this LorentzVector.
        For details, see appendix A.2.2 of Simone Lionetti's PhD thesis.

        :param LorentzVector p: Starting Lorentz vector to define the boost.
        :param LorentzVector q: Target Lorentz vector to define the boost.
        """

        eps = p.eps+q.eps
        p_abs = abs(p)
        q_abs = abs(q)
        assert almost_equal(p.square(), q.square(), rel_tol=eps) or \
               (p.square_almost_zero() and q.square_almost_zero())
        p_vec = p.space()
        q_vec = q.space()
        if almost_equal(p_vec, q_vec, rel_tol=eps):
            return Vector([0 for _ in p_vec])
        n_vec = (q_vec - p_vec).normalize()
        na = LorentzVector([1, ] + list(+n_vec))
        nb = LorentzVector([1, ] + list(-n_vec))
        assert na.square_almost_zero()
        assert nb.square_almost_zero()
        assert almost_equal(na.dot(nb), 2, rel_tol=eps)
        p_plus  = p.dot(nb)
        p_minus = p.dot(na)
        q_plus  = q.dot(nb)
        q_minus = q.dot(na)
        if p_minus/p_abs < eps and q_minus/q_abs < eps:
            if p_plus/p_abs < eps and q_plus/q_abs < eps:
                ratioa = 1
                ratiob = 1
            else:
                ratiob = q_plus / p_plus
                ratioa = 1. / ratiob
        else:
            if p_plus/p_abs < eps and q_plus/q_abs < eps:
                ratioa = q_minus / p_minus
                ratiob = 1. / ratioa
            else:
                ratioa = q_minus / p_minus
                ratiob = q_plus / p_plus
        plus  = self.dot(nb)
        minus = self.dot(na)
        self.__iadd__(((ratiob - 1) * 0.5 * plus ) * na)
        self.__iadd__(((ratioa - 1) * 0.5 * minus) * nb)
        return self

#=========================================================================================
# LorentzVectorDict
#=========================================================================================

class LorentzVectorDict(dict):
    """A simple class wrapping dictionaries that store Lorentz vectors."""

    def to_list(self, ordered_keys=None):
        """Return list copy of self. Notice that the actual values of the keys
        are lost in this process. The user can specify in which order (and potentially which ones)
        the keys must be placed in the list returned."""
        
        if ordered_keys is None:
            return LorentzVectorList(self[k] for k in sorted(self.keys()))
        else:
            return LorentzVectorList(self[k] for k in ordered_keys)

    def to_dict(self):
        """Return a copy of this LorentzVectorDict """

        return LorentzVectorDict(self)

    def to_tuple(self):
        """Return a copy of this LorentzVectorDict as an immutable tuple.
        Notice that the actual values of the keys are lost in this process.
        """

        return tuple( tuple(self[k]) for k in sorted(self.keys()) )

    def __str__(self, n_initial=2):
        """Nice printout of the momenta."""

        # Use padding for minus signs
        def special_float_format(fl):
            return '%s%.16e' % ('' if fl < 0.0 else ' ', fl)

        cols_widths = [4, 25, 25, 25, 25, 25]
        template = ' '.join(
            '%%-%ds' % col_width for col_width in cols_widths
        )
        line = '-' * (sum(cols_widths) + len(cols_widths) - 1)

        out_lines = [template % ('#', ' E', ' p_x', ' p_y', ' p_z', ' M',)]
        out_lines.append(line)
        running_sum = LorentzVector()
        for i in sorted(self.keys()):
            mom = LorentzVector(self[i])
            if i <= n_initial:
                running_sum += mom
            else:
                running_sum -= mom
            out_lines.append(template % tuple(
                ['%d' % i] + [
           special_float_format(el) for el in (list(mom) + [math.sqrt(abs(mom.square()))])
                ]
            ))
        out_lines.append(line)
        out_lines.append(template % tuple(
            ['Sum'] + [special_float_format(el) for el in running_sum] + ['']
        ))

        return '\n'.join(out_lines)

    def boost_to_com(self, initial_leg_numbers):
        """ Boost this kinematic configuration back to its c.o.m. frame given the
        initial leg numbers. This is not meant to be generic and here we *want* to crash
        if we encounter a configuration that is not supposed to ever need boosting in the
        MadNkLO construction.
        """
    
        if len(initial_leg_numbers)==2:
            if __debug__:
                sqrts = math.sqrt((self[initial_leg_numbers[0]]+self[initial_leg_numbers[1]]).square())
                # Assert initial states along the z axis
                assert(abs(self[initial_leg_numbers[0]][1]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[1]][1]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[0]][2]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[1]][2]/sqrts)<1.0e-9)
            # Now send the self back into its c.o.m frame, if necessary
            initial_momenta_summed = self[initial_leg_numbers[0]]+self[initial_leg_numbers[1]]
            sqrts = math.sqrt((initial_momenta_summed).square())
            if abs(initial_momenta_summed[3]/sqrts)>1.0e-9:
                boost_vector = (initial_momenta_summed).boostVector()
                for vec in self.values():
                    vec.boost(-boost_vector)
            if __debug__:
                assert(abs((self[initial_leg_numbers[0]]+self[initial_leg_numbers[1]])[3]/sqrts)<=1.0e-9)
        elif len(initial_leg_numbers)==1:
            if __debug__:
                sqrts = math.sqrt(self[initial_leg_numbers[0]].square())
                assert(abs(self[initial_leg_numbers[0]][1]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[0]][2]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[0]][3]/sqrts)<1.0e-9)
        else:
            raise InvalidOperation('MadNkLO only supports processes with one or two initial states.')

    def get_copy(self):
        """Return a copy that can be freely modified
        without changing the current instance.
        """

        return type(self)((i,LorentzVector(k)) for i,k in self.items())

#=========================================================================================
# LorentzVectorList
#=========================================================================================

class LorentzVectorList(list):
    """A simple class wrapping lists that store Lorentz vectors."""

    def __str__(self, n_initial=2):
        """Nice printout of the momenta."""

        return LorentzVectorDict(
            (i + 1, v) for i, v in enumerate(self)
        ).__str__(n_initial=n_initial)

    def to_list(self):
        """Return list copy of self."""

        return LorentzVectorList(self)

    def to_tuple(self):
        """Return a copy of this LorentzVectorList as an immutable tuple."""

        return tuple( tuple(v) for v in self )

    def to_dict(self):
        """Return a copy of this LorentzVectorList as a LorentzVectorDict."""

        return LorentzVectorDict( (i+1, v) for i, v in enumerate(self) )

    def boost_from_com_to_lab_frame(self, x1, x2, ebeam1, ebeam2):
        """ Boost this kinematic configuration from the center of mass frame to the lab frame
        given specified Bjorken x's x1 and x2.
        This function needs to be cleaned up and built in a smarter way as the boost vector can be written
        down explicitly as a function of x1, x2 and the beam energies.
        """

        if x1 is None: x1 = 1.
        if x2 is None: x2 = 1.

        target_initial_momenta = []
        for i, (x, ebeam) in enumerate(zip([x1, x2],[ebeam1, ebeam2])):
            target_initial_momenta.append(LorentzVector([x*ebeam, 0., 0., math.copysign(x*ebeam, self[i][3])]))
        target_summed = sum(target_initial_momenta)
        source_summed = LorentzVector([2.*math.sqrt(x1*x2*ebeam1*ebeam2),0.,0.,0.])

        # We want to send the source to the target
        for vec in self:
            vec.boost_from_to(source_summed, target_summed)
            #boost_vec = LorentzVector.boost_vector_from_to(source_summed, target_summed)
            #import madgraph.various.misc as misc
            #misc.sprint(boost_vec)
            #vec.boost(boost_vec)

    def boost_to_com(self, initial_leg_numbers):
        """ Boost this kinematic configuration back to its c.o.m. frame given the
        initial leg numbers. This is not meant to be generic and here we *want* to crash
        if we encounter a configuration that is not supposed to ever need boosting in the
        MadNkLO construction.
        """
        # Given that this is a list, we must subtract one to the indices given
        initial_leg_numbers = tuple(n-1 for n in initial_leg_numbers)
        if len(initial_leg_numbers)==2:
            if __debug__:
                sqrts = math.sqrt((self[initial_leg_numbers[0]]+self[initial_leg_numbers[1]]).square())
                # Assert initial states along the z axis
                assert(abs(self[initial_leg_numbers[0]][1]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[1]][1]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[0]][2]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[1]][2]/sqrts)<1.0e-9)
            # Now send the self back into its c.o.m frame, if necessary
            initial_momenta_summed = self[initial_leg_numbers[0]]+self[initial_leg_numbers[1]]
            sqrts = math.sqrt((initial_momenta_summed).square())
            if abs(initial_momenta_summed[3]/sqrts)>1.0e-9:
                boost_vector = (initial_momenta_summed).boostVector()
                for vec in self:
                    vec.boost(-boost_vector)
            if __debug__:
                assert(abs((self[initial_leg_numbers[0]]+self[initial_leg_numbers[1]])[3]/sqrts)<=1.0e-9)
        elif len(initial_leg_numbers)==1:
            if __debug__:
                sqrts = math.sqrt(self[initial_leg_numbers[0]].square())
                assert(abs(self[initial_leg_numbers[0]][1]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[0]][2]/sqrts)<1.0e-9)
                assert(abs(self[initial_leg_numbers[0]][3]/sqrts)<1.0e-9)
        else:
            raise InvalidOperation('MadNkLO only supports processes with one or two initial states.')

    def get_copy(self):
        """Return a copy that can be freely modified
        without changing the current instance.
        """

        return type(self)([LorentzVector(p) for p in self])
