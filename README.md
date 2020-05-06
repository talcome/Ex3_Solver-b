<div dir="rtl" lang="he">

# פותר משוואות - שלב א

האחיין החמוד שלכם בכיתה ה והוא לומד עכשיו לפתור משוואות בנעלם אחד. 
יום אחד הוא בא ושואל אתכם:
"תגידו, אתם יכולים לתכנת את המחשב כך שאני אוכל להגיד לו
'תפתור לי את המשוואה שני איקס פחות ארבע שווה עשר' והוא יגיד לי 'שבע'?"
במטלה זו תוכיחו לו שהתשובה היא "כן". המשימה היא לכתוב מערכת שתאפשר לכם לכתוב:

<div dir='ltr'>

    solve(2*x-4==10);

</div>

ולקבל 7.

כדי לראות איך בדיוק המערכת אמורה לעבוד, ראו בקובץ
[Demo.cpp](Demo.cpp)
המצורף.
שימו לב: המספרים הרשומים בהערה בצד ימין של כל שורה הם התוצאות המספריות הצפויות;
התוצאות לא חייבות להופיע דווקא בפורמט זה
(למשל 7 יכול להופיע גם כ 7.0).
בבדיקות האוטומטיות, רמת הדיוק להשוואת מספרים ממשיים צריכה להיות 0.0001.

בשלב א עליכם לכתוב את הקבצים הדרושים כך שתוכנית הדוגמה תתקמפל.
כמו כן עליכם לכתוב בדיקות-יחידה מקיפות. הפקודות הבאות צריכות לרוץ בלי שגיאות:

<div dir='ltr'>

    make demo && ./demo
	make test && ./test

</div>


עוד פרטים:
* המערכת תומכת בשני סוגי משתנים - ממשי ומרוכב. משוואה במשתנה ממשי מחזירה תוצאה מסוג `double`. משוואה במשתנה מרוכב מחזירה תוצאה מסוג `std::complex`.
* במטלה זו מספיק לפתור משוואות בנעלם אחד, ליניאריות או ריבועיות. אין צורך להתייחס למשוואות בשני נעלמים או יותר, או למשוואות במעלה שלישית או יותר, או לאי-שיוויונים.
במקרה שמישהו מנסה ליצור משוואה כזאת, יש לזרוק חריגה, אבל אין צורך להתייחס לזה בבדיקות האוטומטיות.
* מותר לכתוב מחלקות נוספות או לממש אופרטורים נוספים שאינם נזכרים בתוכנית הראשית, אבל אין צורך לכתוב להם בדיקות.
* אם למשוואה אין פתרון (ממשי או מרוכב לפי סוג המשתנה) יש לזרוק חריגה מתאימה.
* אם למשוואה יש יותר מפתרון אחד, יש להחזיר אחד מהם שרירותית.
* בכל עניין נוסף שאינו מוגדר עד הסוף במטלה, יש לכם חופש לממש לפי שיקול דעתכם, ולכן אין צורך להתייחס לזה בבדיקות. הבדיקות שאתם כותבים צריכות להתמקד בדרישות המטלה בלבד.

</div># solver
