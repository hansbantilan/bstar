import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
from scipy.interpolate import interp1d
import matplotlib
import os
import matplotlib.pyplot as plt
import time


class Complex_Boson_Star:

    edelta_guess = None
    _phi0 = None
    _Dim = None
    _Lambda = None

    verbose = None
    path = None

    _edelta_final = None
    _solution_array = None
    _solution_r_pos = None

    _finished_shooting = False

    def __init__(self, edelta_guess, phi0, Dim, Lambda, verbose=0):

        self.edelta_guess = edelta_guess
        self._phi0 = phi0
        self._Dim = Dim
        self._Lambda = Lambda

        # Will give more messages with increasing value
        self.verbose = verbose

        self.make_file()
        return None

    def print_parameters(self):
        print "----------------------------------------------------"
        print r"The cosmological constant $\Lambda$ ", self._Lambda
        print "The dimension of the problen        ", self._Dim
        print r"Central value of $\phi$             ", self._phi0
        print "----------------------------------------------------"

    def eqns(self, y, r):
        """ Differential equation for scalar fields from arXiv:gr-qc/0309131

        Parameters:
            y (list with reals): current status vector ( a(r), alpha(r), phi(r), pi(r) )
            r (real) : current position

        Returns:
            dydr (list with reals): derivative for y in r

        """
        D = float(self._Dim)
        Lambda = self._Lambda
        edelta, m, phi, pi = y
        # Where edelta  = e^{-\delta}

        F = (1 - 2 * m / r**(D - 3) - 2 * Lambda * r**2 / ((D - 2) * (D - 1)))

        dedeltadr = r * (edelta * pi**2.0 + edelta**(-1) * phi**2 / F**2)
        dmdr = r**(D - 2) * 0.5 * (F * pi**2 + phi **
                                   2 + edelta**(-2) * phi**2 / F)
        dphidr = pi

        dFdr = (-4 * Lambda * r) / ((-2 + D) * (-1 + D)) - 2 * \
            (3 - D) * r**(2 - D) * m - 2 * r**(3 - D) * dmdr

        dpidr = -(phi / (edelta**2 * F**2)) + phi / F - (dedeltadr * \
                  pi) / edelta - (dFdr * pi) / F + (2 * pi) / r - (D * pi) / r
        dydr = [dedeltadr, dmdr, dphidr, dpidr]

        return dydr

    def shoot(self, edelta_at_zero, r, output=False):
        """ Solves differential equation

        Parameters:
            edelta_at_zero (real): The lapse value guess at r = rmin
            r       (real array) : Radial points used for solver
            output  (bool)       : if True outputs whole solution array

        Returns:
            phi_end (real):. The phi value at r = rmax
            or
            sol     (real array) : array containg solution

        """

        # Define initial data vector
        y0 = [edelta_at_zero, 0, self._phi0, 0]
        # Solve differential equaion
        sol = spi.odeint(self.eqns, y0, r)
        phi_end = sol[-1, 2]

        if not output:
            return phi_end
        else:
            return sol

    def radial_walker(self, r_start, r_end, delta_R, N, eps):
        """ Performs shooting for multiple radii rmax shooting process.

        Parameters:
            r_start (real) : first rmax for which shooting is performed
            r_end (real) : maximum rmax for which shooting is performed
            delta_R (real) : stelpsize
            N (real) : number of gridpoints

        Returns:
            alpha0 (real):. alpha0 for rmax
        """
        range_list = np.arange(r_start, r_end, delta_R)
        edelta_guess_tmp = self.edelta_guess

        if self.verbose >= 1:
            print "Shooting started"
        if self.verbose >= 1:
            start = time.time()

        for R_max in range_list:
            r = np.linspace(eps, R_max, N)

            def fun(x): return self.shoot(x, r)
            root = opi.root(fun, edelta_guess_tmp)
            edelta_guess_tmp = root.x

            if self.verbose >= 2:
                print "Edelta at R = eps ", edelta_guess_tmp[0], " with Rmax ", R_max

        if self.verbose >= 1:
            print "Shooting finished in ", time.time() - start, "sec"

        self._finished_shooting = True
        output_solution = True
        self._solution_r_pos = np.linspace(eps, r_end, N)
        self._solution_array = self.shoot(
            edelta_guess_tmp[0],
            self._solution_r_pos,
            output_solution)
        self._edelta_final = edelta_guess_tmp

        return edelta_guess_tmp[0]

    def normalise_edelta(self, sol):
        """ Extractsomega for edelta by the coordinate transformation  t -> omega t

        Parameters:
            sol (real array) : were the sol[:,1] corresponds to edelta^(-1) and
                           and asymtotic value that does not go to 1
        Returns:
            omega (real): frequency of scalar field
            sol (real array) : sol array with fixed edelta
        """

        edelta = 1. / sol[:, 0]
        N = len(edelta)
        omega = edelta[N - 1]
        edelta = edelta / omega
        sol[:, 0] = 1. / edelta
        return omega, sol

    def make_file(self):
        """ Creates Folder for current physics problem if they do not yet exist
        """

        name_Lambda_Dim = "Lambda" + str(self._Lambda) + "D" + str(self._Dim)
        path = name_Lambda_Dim
        if not os.path.exists(path):
            os.mkdir(path)

        name_phi = "phi" + str(self._phi0)
        path = name_Lambda_Dim + "/" + name_phi
        if not os.path.exists(path):
            os.mkdir(path)
            if self.verbose >= 1:
                print "Create Folder with relative ", path
        else:
            if self.verbose >= 1:
                print "Folder with path ", path, " already exists "

        self.path = path

    def get_path(self):
        """ return
              path (string): Realtive path used for outputs
        """
        if self.path is None:
            make_file()
        return self.path

    def get_solution(self):
        """return
             solution_array (real array) : solution array for Rmax
        """
        if self._solution_array is None or self._solution_r_pos is None:
            print("----------------------------------------")
            print("WARNING: SHOOTING HAS NOT BEEN PERFORMED")
            print("----------------------------------------")
            return None
        else:
            return self._solution_r_pos, self._solution_array

    def plot_solution(self):
        """ Prints solution if shooting has been performed already

        """
        if self.path is None:
            make_file()
        if self._solution_array is None or self._solution_r_pos is None:
            print("----------------------------------------")
            print("WARNING: SHOOTING HAS NOT BEEN PERFORMED")
            print("----------------------------------------")
        else:

            if self.verbose >= 1:
                print "Plotting started"
            if self.verbose >= 1:
                start = time.time()

            phi = self._solution_array[:, 2]
            m = self._solution_array[:, 1]
            edelta = 1 / self._solution_array[:, 0]
            r = self._solution_r_pos

            # find 90 % radius of R
            Rguess = 0.01
            maxphi = max(phi)
            phi_tmp_fun = interp1d(r, phi - maxphi * 0.1)
            root = opi.root(phi_tmp_fun, Rguess)
            R90 = root.x[0]

            fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10, 10))
            ax1.plot(r, edelta, 'b', )
            ax2.plot(r, m, 'g')
            ax3.plot(r, phi, 'r')

            ax3.set_xlabel('t')

            ax1.set_ylabel(r'$ e^{\delta (t)}$')
            ax2.set_ylabel('$ m (t)$')
            ax3.set_ylabel(r'$\phi (t)$')

            ax1.set_xlim([0, R90 * 2])
            ax2.set_xlim([0, R90 * 2])
            ax3.set_xlim([0, R90 * 2])

            ax1.grid()
            ax2.grid()
            ax3.grid()

            plt.savefig(self.path + "/overview.png")

            if self.verbose >= 1:
                print "Plotting finished in ", time.time() - start, " sec"
