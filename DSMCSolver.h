#pragma once
#include <cmath>
#include <cstdlib>
#include <vector>
#include <list>
#include <algorithm>
#include <stdexcept>
#include <utility>

#include "DSMCElements.h"
#include "../track3d/ExportOpenFOAM.h"
#include "../track3d/vector3d.hpp"
#include "../track3d/AnsysMesh.h"

namespace EvaporatingParticle
{
using Vector3D = Vector3<double>;

class DSMCSolver : public CAnsysMesh // to use find_elem() and terminate()
{
public:
	DSMCSolver();
	~DSMCSolver();

	double                get_time_step () const;
	DWORD_PTR             get_time_step_ptr() const;
	void                  set_time_step(double time_step_val);

	double                get_exec_time() const;
	DWORD_PTR             get_exec_time_ptr() const;
	void                  set_exec_time(double exec_time_val);

	double                get_molar_mass() const;
	DWORD_PTR             get_molar_mass_ptr() const;
	void                  set_molar_mass(double molar_mass_val);

	double                get_cr_section() const;
	DWORD_PTR             get_cr_section_ptr() const;
	void                  set_cr_section(double cr_section_val);

	double                get_eq_particles() const;
	DWORD_PTR             get_eq_particles_ptr() const;
	void                  set_eq_particles(double eq_particles_val);

	bool		init(CAnsysMesh* mesh, CExportOpenFOAM* exporter);

	bool		init_face(CRegion* pReg, CFace* pFace);


	void		InjectParticles();
	void		MoveParticles();
	void		UpdateCellProperties();
	void		CollideParticles();
	void		SetDataToNodes(CAnsysMesh* mesh);

	std::pair<CRegFacePair, double>		FindIntersection(const CElem3D* elem, CRay& traj, CRegionsCollection& checked_rgns);
	bool								Interact(QuasiParticle* particle, std::pair<CRegFacePair, double> coll_bound);
	void								RegisterParticle(QuasiParticle* particle, const CElem3D* expected_elem, Vector3D expected_pos);
	Vector3D							RandomInletPos(CFace* pFace);

	CElementsCollection*	elements;
	CRegionsCollection*		regions;
	UINT					nStep;

	double			time_step;
	double			exec_time;
	double			molar_mass;
	double			cr_section;
	long double		eq_particles;

	// DEBUG variables
	UINT simCollisions;
};

inline double DSMCSolver::get_time_step() const
{
	return time_step;
}

inline DWORD_PTR DSMCSolver::get_time_step_ptr() const
{
	return (DWORD_PTR)&time_step;
}

inline void DSMCSolver::set_time_step(double time_step_val)
{
	time_step = time_step_val;
}

inline double DSMCSolver::get_exec_time() const
{
	return exec_time;
}

inline DWORD_PTR DSMCSolver::get_exec_time_ptr() const
{
	return (DWORD_PTR)&exec_time;
}

inline void DSMCSolver::set_exec_time(double exec_time_val)
{
	exec_time = exec_time_val;
}

inline double DSMCSolver::get_molar_mass() const
{
	return molar_mass;
}

inline DWORD_PTR DSMCSolver::get_molar_mass_ptr() const
{
	return (DWORD_PTR)&molar_mass;
}

inline void DSMCSolver::set_molar_mass(double molar_mass_val)
{
	molar_mass = molar_mass_val;
}

inline double DSMCSolver::get_cr_section() const
{
	return cr_section;
}

inline DWORD_PTR DSMCSolver::get_cr_section_ptr() const
{
	return (DWORD_PTR)&cr_section;
}

inline void DSMCSolver::set_cr_section(double cr_section_val)
{
	cr_section = cr_section_val;
}

inline double DSMCSolver::get_eq_particles() const
{
	return eq_particles;
}

inline DWORD_PTR DSMCSolver::get_eq_particles_ptr() const
{
	return (DWORD_PTR)&eq_particles;
}

inline void DSMCSolver::set_eq_particles(double eq_particles_val)
{
	cr_section = eq_particles_val;
}

}