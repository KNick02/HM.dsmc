#pragma once

#include <list>
#include <iterator>

#include "../track3d/vector3d.hpp"
#include <../utilities/MemoryPool.h>

namespace EvaporatingParticle
{
using Vector3D = Vector3<double>;

struct dsmcElem;
struct dsmcRegion;
struct QuasiParticle;
struct CElem3D;

struct dsmcElem
{
	std::list<QuasiParticle*>* particles;
	double volume;
	double temp, press;
	double particles_cntr = 0;
	Vector3D average_v;
};

struct dsmcRegion
{
	int nType;
	double fTemp, fPress;
};

struct dsmcFace : public dsmcRegion
{
	const CElem3D* nbrElem = NULL;
};

struct QuasiParticle : public BlockAllocator <QuasiParticle>
{
	Vector3D ExpectedPos(double t_step);

	void Collide(QuasiParticle* qp);

	bool moved;
	Vector3D v, pos;
	const CElem3D* elem;
	std::list<QuasiParticle*>::iterator iter;
};

}