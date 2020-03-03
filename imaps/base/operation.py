"""Base operation."""
import argparse
import inspect


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

    @classmethod
    def cli(cls):
        """Add command line interface."""
        # Create parser
        parser = argparse.ArgumentParser(cls.__doc__)
        for name, parameter in inspect.signature(cls.__init__).parameters.items():
            if name == "self":
                continue
            if parameter.default != inspect._empty:  # pylint: disable=protected-access
                parser.add_argument(f"--{name}", help=name, default=parameter.default)
            else:
                parser.add_argument(name, help=name)

        # Parse arguments
        args = parser.parse_args()

        # Prepare inputs
        inputs = {}
        for name, _ in inspect.signature(cls.__init__).parameters.items():
            if name == "self":
                continue
            inputs[name] = getattr(args, name)

        # Create an instance and run it
        instance = cls(**inputs)
        instance.run()
