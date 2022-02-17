from molSimplify.optimize.coordinates import (Distance, InverseDistance,
                                              Angle, Dihedral, LinearAngle)
import numpy as np
import numdifftools as nd


def test_distances(atol=1e-10):
    e = np.array([1., 2., 3.])
    e /= np.linalg.norm(e)

    d = Distance(0, 1)
    d_inv = InverseDistance(0, 1)
    xyzs = np.zeros((2, 3))
    for r in np.linspace(0.1, 100, 21):
        xyzs[0, :] = 0.5*r*e
        xyzs[1, :] = -0.5*r*e
        assert np.abs(d.value(xyzs) - r) < atol
        assert np.abs(d_inv.value(xyzs) - 1/r) < atol


def test_angle(atol=1e-10):
    a = Angle(1, 0, 2)

    xyzs = np.array([[0., 0., 0.],
                     [1.2, 0., 0.],
                     [0.7, 0., 0.]])
    for theta in [0.01, np.pi/3, np.pi/2, 2*np.pi/3, np.pi - 0.01]:
        xyzs[2, :] = 0.7*np.cos(theta), 0., 0.7*np.sin(theta)
        assert np.abs(a.value(xyzs) - theta) < atol


def test_primitive_derivatives(atol=1e-10):
    xyzs = np.array([[0., 0., 0.],
                     [1.3, 0., 0.],
                     [-0.2, 1.4, 0],
                     [-1.6, 0., 0.],
                     [-0.3, -0.5, 0.9]])

    for i in range(len(xyzs)):
        for j in range(i+1, len(xyzs)):
            r = Distance(i, j)
            dr_ref = nd.Gradient(lambda x: r.value(x.reshape((-1, 3))))
            np.testing.assert_allclose(r.derivative(xyzs), dr_ref(xyzs))

            r_inv = InverseDistance(i, j)
            dr_inv_ref = nd.Gradient(lambda x: r_inv.value(x.reshape((-1, 3))))
            np.testing.assert_allclose(
                r_inv.derivative(xyzs), dr_inv_ref(xyzs), atol=atol)

    a = Angle(1, 0, 2)
    da_ref = nd.Gradient(lambda x: a.value(x.reshape((-1, 3))))
    np.testing.assert_allclose(a.derivative(xyzs), da_ref(xyzs), atol=atol)

    t1 = LinearAngle(1, 0, 3, 0)
    dt1_ref = nd.Gradient(lambda x: t1.value(x.reshape((-1, 3))))
    np.testing.assert_allclose(t1.derivative(xyzs), dt1_ref(xyzs), atol=atol)

    t2 = LinearAngle(1, 0, 3, 1)
    dt2_ref = nd.Gradient(lambda x: t2.value(x.reshape((-1, 3))))
    np.testing.assert_allclose(t2.derivative(xyzs), dt2_ref(xyzs), atol=atol)

    w = Dihedral(1, 0, 2, 4)
    dw_ref = nd.Gradient(lambda x: w.value(x.reshape((-1, 3))))
    np.testing.assert_allclose(w.derivative(xyzs), dw_ref(xyzs), atol=atol)
