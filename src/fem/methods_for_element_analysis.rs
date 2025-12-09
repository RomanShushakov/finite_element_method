use extended_matrix::FloatTrait;

use crate::fem::FEM;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ElementForceComponent {
    ForceR,
    ForceS,
    ForceT,
    MembraneForceR,
    MembraneForceS,
    MembraneForceRS,
    ShearForceRT,
    ShearForceST,
    MomentR,
    MomentS,
    MomentT,
    BendingMomentR,
    BendingMomentS,
    BendingMomentRS,
}

impl<V> FEM<V>
where
    V: FloatTrait<Output = V>,
{
    pub fn extract_elements_analysis_result(
        &self,
    ) -> Result<Vec<(u32, Vec<(ElementForceComponent, V)>)>, String> {
        let mut elements_analysis_result = Vec::new();

        for (truss_element_number, truss_element) in self.get_truss_elements().iter() {
            let element_analysis_result = truss_element.extract_element_analysis_result(
                self.get_nodes(),
                self.get_displacements_vector(),
            )?;
            elements_analysis_result.push((*truss_element_number, element_analysis_result))
        }

        for (beam_element_number, beam_element) in self.get_beam_elements().iter() {
            let element_analysis_result = beam_element.extract_element_analysis_result(
                self.get_nodes(),
                self.get_displacements_vector(),
            )?;
            elements_analysis_result.push((*beam_element_number, element_analysis_result))
        }

        for (plate_element_number, plate_element) in self.get_plate_elements().iter() {
            let element_analysis_result = plate_element.extract_element_analysis_result(
                self.get_nodes(),
                self.get_displacements_vector(),
                self.get_props().get_rel_tol(),
            )?;
            elements_analysis_result.push((*plate_element_number, element_analysis_result))
        }

        Ok(elements_analysis_result)
    }
}
