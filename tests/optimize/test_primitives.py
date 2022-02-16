from molSimplify.optimize.coordinates import (Distance, InverseDistance,
                                              Angle, Dihedral)
import numpy as np
import numdifftools as nd


def test_angle(atol=1e-10):
    a = Angle(1, 0, 2)

    for theta in [0.01, np.pi/3, np.pi/2, 2*np.pi/3, np.pi - 0.01]:
        xyzs = np.array([[0., 0., 0.],
                        [1.2, 0., 0.],
                        [0.7*np.cos(theta), 0., 0.7*np.sin(theta)]])
        assert np.abs(theta - a.q(xyzs)) < atol


def test_primitive_derivatives(atol=1e-10):
    xyzs = np.array([[0., 0., 0.],
                     [1.3, 0., 0.],
                     [-0.2, 1.4, 0],
                     [-1.6, 0., 0.],
                     [-0.3, -0.5, 0.9]])

    for i in range(len(xyzs)):
        for j in range(i+1, len(xyzs)):
            r = Distance(i, j)
            dr_ref = nd.Gradient(lambda x: r.q(x.reshape((-1, 3))))
            np.testing.assert_allclose(r.dq_dx(xyzs), dr_ref(xyzs))

            r_inv = InverseDistance(i, j)
            dr_inv_ref = nd.Gradient(lambda x: r_inv.q(x.reshape((-1, 3))))
            np.testing.assert_allclose(r_inv.dq_dx(xyzs), dr_inv_ref(xyzs))

    a = Angle(1, 0, 2)
    da_ref = nd.Gradient(lambda x: a.q(x.reshape((-1, 3))))
    np.testing.assert_allclose(a.dq_dx(xyzs), da_ref(xyzs), atol=atol)

    w = Dihedral(1, 0, 2, 4)
    dw_ref = nd.Gradient(lambda x: w.q(x.reshape((-1, 3))))
    np.testing.assert_allclose(w.dq_dx(xyzs), dw_ref(xyzs), atol=atol)
