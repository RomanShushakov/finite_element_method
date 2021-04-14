use crate::fem::
    {
        FENode, Truss2n2ip, StiffnessGroup, ElementAnalysisData, Displacements, Beam2n1ipT,
    };

use crate::{ElementsNumbers, ElementsValues};
use crate::extended_matrix::{ExtendedMatrix};

use std::rc::Rc;
use std::cell::RefCell;
use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;

use std::slice::Iter;
use self::FEType::*;


#[derive(Clone, PartialEq, Debug)]
pub enum FEType
{
    Truss2n2ip,
    Beam2n1ipT,
}


impl FEType
{
    pub fn as_str(&self) -> String
    {
        match self
        {
            FEType::Truss2n2ip => String::from("Truss2n2ip"),
            FEType::Beam2n1ipT => String::from("Beam2n1ipT"),
        }
    }


    pub fn iterator() -> Iter<'static, FEType>
    {
        const TYPES: [FEType; 2] =
            [
                Truss2n2ip, Beam2n1ipT,
            ];
        TYPES.iter()
    }
}


pub struct FEData<T, V>
{
    pub number: T,
    pub nodes: Vec<Rc<RefCell<FENode<T, V>>>>,
    pub properties: Vec<V>,
}


impl<T, V> FEData<T, V>
    where T: Debug,
          V: Default + PartialOrd
{
    pub fn check_properties(&self) -> Result<(), String>
    {
        for value in self.properties.iter()
        {
            if *value <= V::default()
            {
                return Err(format!("FEData: All properties values for element {:?} should be \
                    greater than zero!", self.number));
            }
        }
        Ok(())
    }
}


pub trait FiniteElementTrait<T, V>
{
    fn update(&mut self, data: FEData<T, V>) -> Result<(), String>;
    fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>;
    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>;
    fn node_belong_element(&self, node_number: T) -> bool;
    fn refresh(&mut self) -> Result<(), String>;
    fn number_same(&self, number: T) -> bool;
    fn nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool;
    fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>)
        -> Result<ElementAnalysisData<T, V>, String>;
    fn extract_fe_number(&self) -> ElementsNumbers;
    fn extract_nodes_numbers(&self) -> Vec<ElementsNumbers>;
    fn extract_fe_properties(&self) -> Vec<ElementsValues>;
}


struct FECreator<T, V>(T, V);


impl<T, V> FECreator<T, V>
    where T: Copy + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + From<ElementsNumbers> +
             Into<ElementsNumbers> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + PartialOrd +
             Default + Add<Output = T> + 'static,
          V: Copy + From<ElementsValues> + Into<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + 'static,
{
    fn create(fe_type: FEType, data: FEData<T, V>)
        -> Result<Box<dyn FiniteElementTrait<T, V>>, String>
    {
        match fe_type
        {
            FEType::Truss2n2ip =>
                {
                    if data.properties.len() == 3
                    {
                        let truss_element = Truss2n2ip::create(
                            data.number, Rc::clone(&data.nodes[0]),
                            Rc::clone(&data.nodes[1]),
                            data.properties[0], data.properties[1],
                            Some(data.properties[2])
                        )?;
                        Ok(Box::new(truss_element))
                    }
                    else
                    {
                        let truss_element = Truss2n2ip::create(
                            data.number, Rc::clone(&data.nodes[0]),
                            Rc::clone(&data.nodes[1]),
                            data.properties[0], data.properties[1],
                            None
                        )?;
                        Ok(Box::new(truss_element))
                    }
                },
            FEType::Beam2n1ipT =>
                {
                    if data.properties.len() == 3
                    {
                        let beam_element = Beam2n1ipT::create(
                            data.number, Rc::clone(&data.nodes[0]),
                            Rc::clone(&data.nodes[1]),
                            data.properties[0], data.properties[1],
                            Some(data.properties[2])
                        )?;
                        Ok(Box::new(beam_element))
                    }
                    else
                    {
                        let beam_element = Beam2n1ipT::create(
                            data.number, Rc::clone(&data.nodes[0]),
                            Rc::clone(&data.nodes[1]),
                            data.properties[0], data.properties[1],
                            None
                        )?;
                        Ok(Box::new(beam_element))
                    }
                }
        }
    }
}


pub struct FiniteElement<T, V>
{
    element_type: FEType,
    element: Box<dyn FiniteElementTrait<T, V>>,
}


impl<T, V> FiniteElement<T, V>
    where T: Copy + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + From<ElementsNumbers> +
             Into<ElementsNumbers> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + PartialOrd +
             Default + Add<Output = T> + 'static,
          V: Copy + From<ElementsValues> + Into<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + 'static,
{
    pub fn create(fe_type: FEType, data: FEData<T, V>) -> Result<Self, String>
    {
        let element = FECreator::create(fe_type.clone(), data)?;
        Ok(FiniteElement { element_type: fe_type, element })
    }


    pub fn update(&mut self, data: FEData<T, V>) -> Result<(), String>
    {
        self.element.update(data)?;
        Ok(())
    }


    pub fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>
    {
        let stiffness_matrix = self.element.extract_stiffness_matrix()?;
        Ok(stiffness_matrix)
    }


    pub fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>
    {
        self.element.extract_stiffness_groups()
    }


    pub fn node_belong_element(&self, node_number: T) -> bool
    {
        self.element.node_belong_element(node_number)
    }


    pub fn refresh(&mut self) -> Result<(), String>
    {
        self.element.refresh()?;
        Ok(())
    }


    pub fn type_same(&self, element_type: &FEType) -> bool
    {
        self.element_type == *element_type
    }


    pub fn number_same(&self, number: T) -> bool
    {
        self.element.number_same(number)
    }


    pub fn nodes_numbers_same(&self, nodes_numbers: Vec<T>) -> bool
    {
        self.element.nodes_numbers_same(nodes_numbers)
    }


    pub fn extract_element_analysis_data(&self, global_displacements: &Displacements<T, V>)
        -> Result<ElementAnalysisData<T, V>, String>
    {
        let element_analysis_data =
            self.element.extract_element_analysis_data(global_displacements)?;
        Ok(element_analysis_data)
    }


    pub fn extract_fe_type(&self) -> FEType
    {
        self.element_type.clone()
    }


    pub fn extract_fe_number(&self) -> ElementsNumbers
    {
        self.element.extract_fe_number()
    }


    pub fn extract_nodes_numbers(&self) -> Vec<ElementsNumbers>
    {
        self.element.extract_nodes_numbers()
    }


    pub fn extract_fe_properties(&self) -> Vec<ElementsValues>
    {
        self.element.extract_fe_properties()
    }
}
