from DescriptorElement import DescriptorElement

__author__ = 'modlab'

# One could do the property access with lambda function and give the charge as default lambda or s.th. This would easily
# Allow for application of inheritance


class ChargeShell(DescriptorElement):
    description = "Vector of charge per radial charge shells. Shells are defined per bonds to starting atom."

    def __init__(self, **options):
        self.n_shells = options.pop('n_shells')
        self.charge_type = options.pop('charge_type')
        self.shell_differences = options.pop('shell_differences', False)
        super(ChargeShell, self).__init__(**options)

    def calculate_elements(self, atom):
        charge = self.charge_type
        last_shell = []
        tmp_shell = []
        descriptor_elements = []

        contained_atoms = set()
        contained_atoms.add(atom)
        last_shell.append(atom)
        # set the first element of descriptor elements to starting atom charge
        last_shell_charge = float(atom.GetProp(charge))
        descriptor_elements.append(last_shell_charge)

        for shell in range(self.n_shells):
            for ls_atom in last_shell:
                for neighbour_atom in ls_atom.GetNeighbors():
                    if neighbour_atom.GetAtomicNum() != 1:
                        tmp_shell.append(neighbour_atom)
            current_shell = [atm for atm in tmp_shell if atm not in contained_atoms]
            contained_atoms.update(current_shell)
            # Append the partial charges to the list that is later added to the descriptor vector
            shell_charges = [float(tmp_atm.GetProp(charge)) for tmp_atm in current_shell]
            current_shell_charge = sum(shell_charges) / len(shell_charges)
            if self.shell_differences:
                descriptor_elements += (current_shell_charge - last_shell_charge)
            else:
                descriptor_elements += [current_shell_charge]
            tmp_shell = []
            last_shell = current_shell
            last_shell_charge = current_shell_charge

        assert len(descriptor_elements) == self.n_shells + 1
        return descriptor_elements


class MassShell(DescriptorElement):
    description = "Vector of Mass per radial shell. Gives mass and then should somehow account for the steric bulk"

    def __init__(self, **options):
        self.n_shells = options.pop('n_shells')
        self.shell_differences = options.pop('shell_differences', False)
        super(MassShell, self).__init__(**options)

    def calculate_elements(self, atom):
        last_shell = []
        tmp_shell = []
        descriptor_elements = []

        contained_atoms = set()
        contained_atoms.add(atom)
        last_shell.append(atom)
        # set the first element of descriptor elements to starting atom charge
        last_shell_mass = float(atom.GetMass())
        descriptor_elements.append(last_shell_mass)

        for shell in range(self.n_shells):
            for ls_atom in last_shell:
                for neighbour_atom in ls_atom.GetNeighbors():
                    if neighbour_atom.GetAtomicNum() != 1:
                        tmp_shell.append(neighbour_atom)
            current_shell = [atm for atm in tmp_shell if atm not in contained_atoms]
            contained_atoms.update(current_shell)
            # Append the partial charges to the list that is later added to the descriptor vector
            shell_masses = [float(tmp_atm.GetMass()) for tmp_atm in current_shell]
            # The mass of the shell is summed and NOT divided by the atom numbers
            current_shell_mass = sum(shell_masses)
            if self.shell_differences:
                descriptor_elements += (current_shell_mass - last_shell_mass)
            else:
                descriptor_elements += [current_shell_mass]
            tmp_shell = []
            last_shell = current_shell
            last_shell_mass = current_shell_mass

        assert len(descriptor_elements) == self.n_shells + 1
        return descriptor_elements
