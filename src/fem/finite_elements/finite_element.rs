use crate::fem::{FeNode, Truss2n2ip};
use crate::{ElementsNumbers, ElementsValues};
use crate::extended_matrix::{ExtendedMatrix, MatrixElementPosition};

use std::rc::Rc;
use std::cell::RefCell;
use std::ops::{Sub, Div, Rem, SubAssign, Mul, Add, AddAssign, MulAssign};
use std::hash::Hash;
use std::fmt::Debug;


#[derive(Debug)]
pub enum StiffnessType
{
    Kuu,
    Kuth,
    Kthu,
    Kthth,
}


#[derive(Debug)]
pub struct StiffnessGroup<T>
{
    pub stiffness_type: StiffnessType,
    pub number_1: T,
    pub number_2: T,
    pub positions: Vec<MatrixElementPosition<T>>
}


#[derive(Clone)]
pub enum FEType
{
    Truss2n2ip
}


pub struct FEData<T, V>
{
    pub number: T,
    pub nodes: Vec<Rc<RefCell<FeNode<T, V>>>>,
    pub properties: Vec<V>,
}


pub trait FiniteElementTrait<T, V>
{
    fn update(&mut self, data: FEData<T, V>) -> Result<(), &str>;
    fn extract_stiffness_matrix(&self) -> Result<ExtendedMatrix<T, V>, &str>;
    fn extract_stiffness_groups(&self) -> Vec<StiffnessGroup<T>>;
    fn node_belong_element(&self, node_number: T) -> bool;
    fn refresh(&mut self) -> Result<(), &str>;
    fn number_same(&self, number: T) -> bool;
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
    fn create<'a>(fe_type: FEType, data: FEData<T, V>)
        -> Result<Box<dyn FiniteElementTrait<T, V>>, &'a str>
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
                }
        }
    }
}


pub struct FiniteElement<T, V>
{
    pub element_type: FEType,
    pub element: Box<dyn FiniteElementTrait<T, V>>,
}


impl<T, V> FiniteElement<T, V>
    where T: Copy + Sub<Output = T> + Div<Output = T> + Rem<Output = T> + From<ElementsNumbers> +
             Into<ElementsNumbers> + Eq + Hash + SubAssign + Debug + Mul<Output = T> + PartialOrd +
             Default + Add<Output = T> + 'static,
          V: Copy + From<ElementsValues> + Into<ElementsValues> + Sub<Output = V> + Default +
             Mul<Output = V> + Add<Output = V> + Div<Output = V> + PartialEq + Debug + AddAssign +
             MulAssign + SubAssign + 'static,
{
    pub fn create<'a>(fe_type: FEType, data: FEData<T, V>) -> Result<Self, &'a str>
    {
        let element_type = fe_type.clone();
        let element = FECreator::create(fe_type, data)?;
        Ok(FiniteElement { element_type, element })
    }


    pub fn update(&mut self, data: FEData<T, V>) -> Result<(), &str>
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


    pub fn refresh(&mut self) -> Result<(), &str>
    {
        self.element.refresh()?;
        Ok(())
    }


    pub fn number_same(&self, number: T) -> bool
    {
        self.element.number_same(number)
    }
}