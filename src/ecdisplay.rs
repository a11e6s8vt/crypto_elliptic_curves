pub trait EcDisplay {
    fn format_curve(&self) -> anyhow::Result<String>;
}
