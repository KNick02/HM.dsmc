#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <list>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>

#include "stdafx.h"
#include "DSMCSolver.h"
#include "DSMCElements.h"
#include "../track3d/vector3d.hpp"
#include "../track3d/ExportOpenFOAM.h"

#define		R		8.314
#define		PI		3.1416

using std::sqrt;
using std::cos;
using std::sin;

namespace EvaporatingParticle
{
	DSMCSolver::DSMCSolver()
	{
		
	}

	DSMCSolver::~DSMCSolver()
	{
		for (size_t i = 0; i < elements->size(); i++)
		{
			std::list<QuasiParticle*>* ptrList = elements->at(i)->particles;
			for (auto ptrParticle : *ptrList)
				delete ptrParticle;
			delete ptrList;
		}
	}

	bool DSMCSolver::init(CAnsysMesh* mesh, CExportOpenFOAM* exporter)
	{
		elements	= &(mesh->get_elems());
		regions		= &(mesh->get_regions(false));

		size_t conditions_count = exporter->get_bc_count();
		if (!conditions_count)
			return false;

		for (size_t i = 0; i < conditions_count; i++)
		{
			CBoundaryConditions* cond = exporter->get_bound_cond(i);
			for (const std::string& regName : cond->vRegNames)
			{
				CRegion* pReg = get_region(regName);

				pReg->nType		= cond->nType;
				pReg->fPress	= cond->fPress;
				pReg->fTemp		= cond->fTemp;

				for (CFace* pFace : pReg->vFaces)
					if (!init_face(pReg, pFace))
						return false;
			}
		}

		for (size_t i = 0; i < elements->size(); i++)
		{
			elements->at(i)->particles = new std::list<QuasiParticle*>;
			CBox& b = elements->at(i)->box;				// determine the volume of Element, to be corrected
			elements->at(i)->volume = (b.vMax.x - b.vMin.x) 
				* (b.vMax.y - b.vMin.y)
				* (b.vMax.z - b.vMin.z) / 6;
		}

		return true;
	}

	bool DSMCSolver::init_face(CRegion* pReg, CFace* pFace)
	{
		pFace->fPress	= pReg->fPress;
		pFace->fTemp	= pReg->fTemp;
		pFace->nType	= pReg->nType;

		// search of the neighbouring CElem3D
		unsigned match_nodes;
		CElementsCollection p0_elems = pFace->p0->get_nbr_elems();

		for (const CElem3D* pElem : p0_elems)
		{
			match_nodes = 0;
			for (CNode3D* pNode : pElem->get_nodes())
				if ((pNode == pFace->p0) || (pNode == pFace->p1) || (pNode == pFace->p2))
					match_nodes++;

			if (match_nodes == 3)
			{
				pFace->nbrElem = pElem;
				return true;
			}
		}

		return false;		// no neighbour Elem found
	}



	void DSMCSolver::InjectParticles()
	{
		std::string sFileName = "logInjectParticles.data";
		FILE* pStream;
		errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));

		for (CRegion* pReg : *regions)
		{
			if (pReg->nType == bcWall)
				continue;

			for (CFace* pFace : pReg->vFaces)
			{
				const CElem3D* pElem = pFace->nbrElem;
				double inflow = (pFace->fPress / (R * pFace->fTemp)) * sqrt(8 * R * pFace->fTemp / molar_mass / PI) * pFace->square() * time_step / 4;
				UINT nParticles = inflow / eq_particles;
				fprintf(pStream, "%s : %lg \n", pReg->sName, inflow);
				std::default_random_engine gen;
				std::gamma_distribution<> d(1, 2);

				Vector3D et1 = pFace->p1->pos - pFace->p0->pos;
				Vector3D et2 = et1 * pFace->norm;

				for (size_t j = 0; j < nParticles; j++)
				{
					QuasiParticle* new_particle = new QuasiParticle();

					new_particle->iter	= pElem->particles->insert(pElem->particles->end(), new_particle);
					new_particle->moved = false;
					new_particle->pos	= RandomInletPos(pFace);
					new_particle->elem	= pElem;

					double rg1 = d(gen), rg2 = d(gen);
					double r = (double)(rand()) / (double)(RAND_MAX);

					new_particle->v		= sqrt(R * pFace->fTemp / molar_mass) * 
						(pFace->norm * sqrt(rg1) + (et1 * cos(2 * PI * r) + et2 * sin(2 * PI * r)) * sqrt(rg2));
				}
			}
		}
		fclose(pStream);
	}


	void DSMCSolver::MoveParticles()
	{
		// should add division into threads here
		for (const CElem3D* pElem : *elements)
			for (QuasiParticle* particle : *(pElem->particles))
			{
				if (particle->moved)			// skip just added particles
				{
					particle->moved = false;
					continue;
				}

				Vector3D expected_pos = particle->ExpectedPos(time_step);
				const CElem3D* expected_elem = find_elem(pElem, expected_pos);

				if (expected_elem == NULL)		// particle interacts with a bound
				{
					CRay traj(particle->pos, expected_pos);
					const CElem3D* prev_cell = pElem;
					CRegionsCollection checked_rgns;

					std::pair<CRegFacePair, double> coll_bound = FindIntersection(prev_cell, traj, checked_rgns);
					bool collision = Interact(particle, coll_bound);
				}

				else
					RegisterParticle(particle, expected_elem, expected_pos);
			}
	}

	void DSMCSolver::UpdateCellProperties()
	{
		std::string sFileName = "logUpdateCellProperties.data";
		FILE* pStream;
		errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));

		nStep++;

		for (CElem3D* pElem : *elements)
		{
			Vector3D elemV(0, 0, 0);
			double elemE = 0, elemP = 0;
			pElem->particles_cntr += pElem->particles->size();

			for (QuasiParticle* qp : *(pElem->particles))
			{
				elemE += qp->v & qp->v;
				elemV += qp->v;
			}

			pElem->average_v		= elemV / pElem->particles->size();
			pElem->temp				= elemE * molar_mass / (3 * R * pElem->particles->size());
			pElem->press			= (pElem->particles->size() * eq_particles * R * pElem->temp) / pElem->volume;
			fprintf(pStream, "%u %f \n", pElem->particles->size(), pElem->temp);
		}
		fclose(pStream);
	}

	void DSMCSolver::CollideParticles()
	{
		std::string sFileName = "logCollideParticles.data";
		FILE* pStream;
		errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));

		for (CElem3D* pElem : *elements)
		{
			std::list<QuasiParticle*>& pList = *(pElem->particles);

			UINT nParticles = pList.size();
			UINT nCollisions = 0.5 * nParticles * (pElem->particles_cntr / nStep) * eq_particles *
				1e-20 * cr_section * sqrt(2 * R * pElem->temp / molar_mass) * time_step / pElem->volume;
		
			if (nCollisions > 0.5 * nParticles * (nParticles - 1))
				throw std::runtime_error("DSMC::CollideParticles(): too many collisions");

			fprintf(pStream, "%u ", nCollisions);
			std::set<std::set<UINT>> pairs;
			for (UINT i = 0; i < nCollisions;)
			{
				UINT f_ind = rand() % nParticles, s_ind = rand() % nParticles;
				if (f_ind == s_ind)
					continue;

				i += pairs.insert({ f_ind, s_ind }).second;
			}

			for (const std::set<UINT>& pair : pairs)
			{
				auto pair_it = pair.begin();
				UINT p1_ind = *pair_it, p2_ind = *(++pair_it);

				auto pList_it = pList.begin();
				std::advance(pList_it, p1_ind);
				QuasiParticle* p1 = *pList_it;

				pList_it = pList.begin();
				std::advance(pList_it, p2_ind);
				QuasiParticle* p2 = *pList_it;

				p1->Collide(p2);
			}
		}
		fclose(pStream);
	}

	void DSMCSolver::SetDataToNodes(CAnsysMesh* mesh)
	{
		set_job_name("Setting data to mesh nodes...");

		CNodesVector& nodes = mesh->get_nodes();
		for (size_t i = 0; i < nodes.size(); i++)
		{
			set_progress(int(100. * i / nodes.size()));

			CNode3D& node = nodes.at(i);
			if (node.vNbrElems.size() == 0)
				continue;

			double nodeTemp = 0;
			double nodePress = 0;

			for (UINT elemInd : node.vNbrElems)
			{
				nodeTemp += elements->at(elemInd)->temp;
				nodePress += elements->at(elemInd)->press;
			}
			
			nodeTemp /= node.vNbrElems.size();
			nodePress /= node.vNbrElems.size();

			node.temp = nodeTemp;
			node.press = nodePress;
		}
	}

	std::pair<CRegFacePair, double> DSMCSolver::FindIntersection(const CElem3D* elem, CRay& traj, CRegionsCollection& checked_rgns)
	{
		CNodesCollection nbr_nodes = elem->get_nodes();
		double dist;
		UINT face_id;

		for (size_t node = 0; node < nbr_nodes.size(); node++)
		{
			CFaceIndices bnd_faces = nbr_nodes.at(node)->vNbrFaces;
			for (size_t face = 0; face < bnd_faces.size(); face++)
			{
				if (std::find(checked_rgns.begin(), checked_rgns.end(), regions->at(bnd_faces.at(face).nReg))
					== checked_rgns.end() && regions->at(bnd_faces.at(face).nReg)->intersect(traj, dist, face_id))
					return std::make_pair(CRegFacePair(bnd_faces.at(face).nReg, face_id), dist);
				checked_rgns.push_back(regions->at(bnd_faces.at(face).nReg));
			}
		}

		for (size_t region = 0; region < regions->size(); region++)
		{
			if (std::find(checked_rgns.begin(), checked_rgns.end(), regions->at(region))
				== checked_rgns.end() && regions->at(region)->intersect(traj, dist, face_id))
				return std::make_pair(CRegFacePair(region, face_id), dist);
			checked_rgns.push_back(regions->at(region));
		}

		throw std::runtime_error("DSMCSolver::FindIntersection(): no intersection found");
	}

	bool DSMCSolver::Interact(QuasiParticle* particle, std::pair<CRegFacePair, double> coll_bound)
	{
		if (regions->at(coll_bound.first.nReg)->nType != 0)		// flies out of system
		{
			particle->elem->particles->erase(particle->iter);
			delete particle;
			return false;
		}

		double time_to_collide = coll_bound.second / particle->v.length();
		Vector3D n = regions->at(coll_bound.first.nReg)->vFaces.at(coll_bound.first.nFace)->norm;

		particle->pos += particle->v * time_to_collide;
		particle->v -= 2 * (particle->v & n) * n;		// specular reflection

		Vector3D expected_pos = particle->pos + particle->v * (time_step - time_to_collide);
		const CElem3D* expected_elem = find_elem(particle->elem, expected_pos);

		if (expected_elem == NULL)		// double reflection probably, do not compute the second one
		{
			expected_pos = particle->pos;
			expected_elem = find_elem(particle->elem, particle->pos);
			if (expected_elem == NULL)
				throw std::runtime_error("DSMCSolver::Interact(): particle is lost after rebound");
		}

		RegisterParticle(particle, expected_elem, expected_pos);
		return true;
	}

	void DSMCSolver::RegisterParticle(QuasiParticle* particle, const CElem3D* expected_elem, Vector3D expected_pos)
	{
		particle->pos = expected_pos;
		if (expected_elem != particle->elem)
		{
			expected_elem->particles->splice(
				expected_elem->particles->end(), *(particle->elem->particles), particle->iter);
			particle->iter = expected_elem->particles->end();
		}
		particle->elem = expected_elem;
		particle->moved = true;
	}

	Vector3D DSMCSolver::RandomInletPos(CFace* pFace)
	{
		Vector3F a = pFace->p1->pos - pFace->p0->pos;
		Vector3F b = pFace->p2->pos - pFace->p0->pos;

		float k1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float k2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

		if (k1 + k2 > 1)
		{
			k1 = 1 - k1;
			k2 = 1 - k2;
		}

		return k1 * a + k2 * b;
	}
}