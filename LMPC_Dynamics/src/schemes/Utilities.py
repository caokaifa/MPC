import numpy as np
import pdb


def Regression(x, u, lamb):
    Y = x[2:x.shape[0], :]
    X = np.hstack((x[1:(x.shape[0] - 1), :], u[1:(x.shape[0] - 1), :]))

    Q = np.linalg.inv(np.dot(X.T, X) + lamb * np.eye(X.shape[1]))
    b = np.dot(X.T, Y)
    W = np.dot(Q, b)

    A = W.T[:, 0:6]
    B = W.T[:, 6:8]

    ErrorMatrix = np.dot(X, W) - Y
    ErrorMax = np.max(ErrorMatrix, axis=0)
    ErrorMin = np.min(ErrorMatrix, axis=0)
    Error = np.vstack((ErrorMax, ErrorMin))

    return A, B, Error


def Curvature(s, PointAndTangent):
    TrackLength = PointAndTangent[-1,3]+PointAndTangent[-1,4]

    while s > TrackLength:
        s = s - TrackLength

    index = np.all([[s >= PointAndTangent[:, 3]], [s < PointAndTangent[:, 3] + PointAndTangent[:, 4]]], axis=0)

    i = int(np.where(np.squeeze(index))[0])
    curvature = PointAndTangent[i, 5]

    return curvature


def getAngle(s, epsi, PointAndTangent):
    TrackLength = PointAndTangent[-1, 3]+PointAndTangent[-1, 4]

    while s > TrackLength:
        s = s - TrackLength

    index = np.all([[s >= PointAndTangent[:, 3]], [s < PointAndTangent[:, 3] + PointAndTangent[:, 4]]], axis=0)
    i = int(np.where(np.squeeze(index))[0])

    if i > 0:
        ang = PointAndTangent[i - 1, 2]
    else:
        ang = 0

    if PointAndTangent[i, 5] == 0:
        r= 0
    else:
        r = 1 / PointAndTangent[i, 5]

    if r == 0:
        angle_at_s = ang + epsi
    else:
        cumulative_s = PointAndTangent[i, 3]
        relative_s = s - cumulative_s
        spanAng = relative_s / np.abs(r)
        psi = wrap(ang + spanAng * np.sign(r))
        # pdb.set_trace()
        angle_at_s = psi + epsi

    return angle_at_s


def wrap(angle):
    if angle < -np.pi:
        w_angle = 2 * np.pi + angle
    elif angle > np.pi:
        w_angle = angle - 2 * np.pi
    else:
        w_angle = angle

    return w_angle
