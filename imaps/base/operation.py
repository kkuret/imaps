"""Base operation."""


class BaseOperation:
    """Base operation."""

    def validate_inputs(self):
        """Validate inputs."""
        raise NotImplementedError("Overwrite this in subclasses.")

    def main(self):
        """Run."""
        raise NotImplementedError("Overwrite this in subclasses.")

    def run(self):
        """Run."""
        self.validate_inputs()
        self.main()
