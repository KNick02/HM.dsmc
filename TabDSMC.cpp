
#include "stdafx.h"

#include "../ui/PropertiesWnd.h"
#include "../ui/ResponseProperty.h"
#include "../ui/Button.h"
#include "../ui/ParticleTracking.h"
#include "DSMCSolver.h"


void CPropertiesWnd::add_dsmc_ctrls()
{
	EvaporatingParticle::DSMCSolver* pSolver = CParticleTrackingApp::Get()->GetSolver();

	CMFCPropertyGridProperty* pTimeGroup = new CMFCPropertyGridProperty(_T("DSMC computation parameters"));

	CMFCPropertyGridProperty* pTStepProp = new CMFCPropertyGridProperty(_T("Time Step, s"), COleVariant(pSolver->get_time_step()), _T("Defines the time step of the set of three simulation stages"), pSolver->get_time_step_ptr());
	pTimeGroup->AddSubItem(pTStepProp);

	CMFCPropertyGridProperty* pExecTimeProp = new CMFCPropertyGridProperty(_T("Execution time, s"), COleVariant(pSolver->get_exec_time()), _T("Defines the time of general DSMC simulation"), pSolver->get_exec_time_ptr());
	pTimeGroup->AddSubItem(pExecTimeProp);

	m_wndPropList.AddProperty(pTimeGroup);

	CMFCPropertyGridProperty* pParamGroup = new CMFCPropertyGridProperty(_T("Particle parameters"));

	CMFCPropertyGridProperty* pMassProp = new CMFCPropertyGridProperty(_T("Molar mass, kg/mol"), COleVariant(pSolver->get_molar_mass()), _T("Value of considered particles' molar mass"), pSolver->get_molar_mass_ptr());
	pParamGroup->AddSubItem(pMassProp);

	CMFCPropertyGridProperty* pCrossSecProp = new CMFCPropertyGridProperty(_T("Cross section, A * A"), COleVariant(pSolver->get_cr_section()), _T("Value of considered particles' molar mass"), pSolver->get_cr_section_ptr());
	pParamGroup->AddSubItem(pCrossSecProp);

	CMFCPropertyGridProperty* pEqPartProp = new CMFCPropertyGridProperty(_T("Number of equivalent molecules"), COleVariant(pSolver->get_eq_particles()), _T("Defines the number of real gas molecules which represents one quasi-particle"), pSolver->get_eq_particles_ptr());
	pParamGroup->AddSubItem(pEqPartProp);

	m_wndPropList.AddProperty(pParamGroup);
}


void CPropertiesWnd::set_dsmc_data()
{
	EvaporatingParticle::DSMCSolver* pSolver = CParticleTrackingApp::Get()->GetSolver();

	CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pSolver->get_time_step_ptr());
	if (pProp != NULL)
		pSolver->set_time_step(pProp->GetValue().dblVal);

	pProp = m_wndPropList.FindItemByData(pSolver->get_exec_time_ptr());
	if (pProp != NULL)
		pSolver->set_exec_time(pProp->GetValue().dblVal);

	pProp = m_wndPropList.FindItemByData(pSolver->get_molar_mass_ptr());
	if (pProp != NULL)
		pSolver->set_molar_mass(pProp->GetValue().dblVal);

	pProp = m_wndPropList.FindItemByData(pSolver->get_cr_section_ptr());
	if (pProp != NULL)
		pSolver->set_cr_section(pProp->GetValue().dblVal);

	pProp = m_wndPropList.FindItemByData(pSolver->get_eq_particles_ptr());
	if (pProp != NULL)
		pSolver->set_eq_particles(pProp->GetValue().dblVal);
}


void CPropertiesWnd::update_dsmc_ctrls()
{
	EvaporatingParticle::DSMCSolver* pSolver = CParticleTrackingApp::Get()->GetSolver();
	bool bEnable = !m_bBusy;

	CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pSolver->get_time_step_ptr());
	if (pProp != NULL)
		pProp->Enable(bEnable);

	pProp = m_wndPropList.FindItemByData(pSolver->get_exec_time_ptr());
	if (pProp != NULL)
		pProp->Enable(bEnable);

	pProp = m_wndPropList.FindItemByData(pSolver->get_molar_mass_ptr());
	if (pProp != NULL)
		pProp->Enable(bEnable);

	pProp = m_wndPropList.FindItemByData(pSolver->get_cr_section_ptr());
	if (pProp != NULL)
		pProp->Enable(bEnable);
}