/** 
 * JanetTest.java is an application that compiles and runs
 * under J2SE 5.0. It requires no other files.
 */

public class JanetTest {
    public static void main(String[] args) {

        //a number
        float number1 = 37;
        float number2 = 0.00000001f;
        
        System.out.println("Variable values...");
        System.out.println("    number1 = " + number1);
        System.out.println("    number2 = " + number2);

        System.out.println("tests...");
       
        System.out.print(String.valueOf(number2));

        System.out.println("");

       
        if ((number1 - 0.00003)< 0) {
            System.out.println("number1 is close to zero");
        } else {
            System.out.println("number1 is not close to zero");
        }

        if (  (number2 - 0.00003) < 0) {
            System.out.println("number2 is close to zero");
        } else {
            System.out.println("number2 is not close to zero");
        }

    }
}
