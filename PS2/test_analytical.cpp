#include "analytical_integration.h"

int main() {

    cout << "Analytical Integration" << endl;
    // Create shells for testing
    Shell s_shell_origin = Shell(vec({0, 0, 0}), 1.0, 0);
    Shell s_shell_offset = Shell(vec({1, 1, 1}), 1.0, 0);
    Shell p_shell_origin = Shell(vec({0, 0, 0}), 1.0, 1);
    Shell p_shell_offset = Shell(vec({1, 1, 1}), 1.0, 1);

    // Create shell pairs for testing
    ShellPair ss_origin(s_shell_origin, s_shell_origin);
    ShellPair ss_offset(s_shell_origin, s_shell_offset);
    ShellPair sp_origin(s_shell_origin, p_shell_origin);
    ShellPair sp_offset(s_shell_origin, p_shell_offset);
    ShellPair pp_origin(p_shell_origin, p_shell_origin);
    ShellPair pp_offset(p_shell_origin, p_shell_offset);

    // Test 1D overlap integrals
    cout << "Testing 1D Overlap Integrals" << endl;
    cout << "----------------------------" << endl;
    cout << "- SS Integral at Origin: " << setw(20) << ss_origin.overlapIntegral1D(0, 0, 0) << endl;
    cout << "- SP Integral at Origin: " << setw(20) << sp_origin.overlapIntegral1D(0, 0, 1) << endl;
    cout << "- SS Integral with Offset of 1: " << setw(13) << ss_offset.overlapIntegral1D(0, 0, 0) << endl;
    cout << "- SP Integral with Offset of 1: " << setw(13) << sp_offset.overlapIntegral1D(0, 0, 1) << endl;
    cout << endl;

    // Test 3D overlap integrals
    cout << "Testing 3D Overlap Integrals" << endl;
    cout << "----------------------------" << endl;
    cout << "SS Integral at Origin: " << endl;
    cout << ss_origin.overlapIntegral3D() << endl;
    cout << "SP Integral at Origin: " << endl;
    cout << sp_origin.overlapIntegral3D() << endl;
    cout << "PP Integral at Origin: " << endl;
    cout << pp_origin.overlapIntegral3D() << endl;
    cout << "SS Integral with Offset of 1: " << endl;
    cout << ss_offset.overlapIntegral3D() << endl;
    cout << "SP Integral with Offset of 1: " << endl;
    cout << sp_offset.overlapIntegral3D() << endl;
    cout << "PP Integral with Offset of 1: " << endl;
    cout << pp_offset.overlapIntegral3D() << endl;

    return 0;
}