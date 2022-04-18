#include <list>

#include "stdafx.h"
#include "DSMCElements.h"
#include "../track3d/vector3d.hpp"

namespace EvaporatingParticle
{
	Vector3D QuasiParticle::ExpectedPos(double t_step)
	{
		return pos + v * t_step;
	}

	void QuasiParticle::Collide(QuasiParticle* qp)
	{
		double alpha_1 = ((double)std::rand() / (double)RAND_MAX);
		double alpha_2 = ((double)std::rand() / (double)RAND_MAX);
		double cos_th = 1 - 2 * alpha_1, sin_th = std::sqrt(1 - cos_th * cos_th);

		Vector3D norm(cos_th, sin_th * std::cos(2 * 3.1416 * alpha_2), sin_th * std::sin(2 * 3.1416 * alpha_2));
		v += norm * ((qp->v - v) & norm);
		qp->v -= norm * ((qp->v - v) & norm);
	}

}